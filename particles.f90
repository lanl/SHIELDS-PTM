module particles

use global
use rksuite_90
use fields
use interpolation

integer, parameter, private :: irand = 2001

integer :: idist, idens
real(dp) :: x0, y0, z0, E0, alpha, vtperp, vtpara, vz0, vp0
real(dp) :: xsMin, xsMax, ysMin, ysMax, zsMin, zsMax, mltMin, mltMax, r0

! Flux mapping
integer :: nEnergies, nPitchAngles
real(dp) :: EkevMin, EkevMax, PitchAngleDegMin, PitchAngleDegMax
real(dp), dimension(:), allocatable :: pitchAngles, energies

contains

  subroutine particles_initialize()

  implicit none

  integer :: nseed, lun, ierr
  real(dp) :: gam, vt
  integer :: i

  ! Initialize the system RNG
  call random_seed(size=nseed)
  call random_seed(put=(/(irand+i,i=0,nseed-1)/))

  lun = get_open_unit()
  open(unit=lun,file='ptm_input/dist_velocity_'//id_string//'.txt',action='read',status='old',iostat=ierr)
  call assert(ierr==0,'particles_initialize','error opening dist_velocity_'//id_string//'.txt')

  read(lun,*) idist

  select case(idist)
    case(1) ! Beam
      read(lun,*) E0
      read(lun,*) alpha
      gam = 1.0+E0/(mass*mc2)
      vt = ckm*sqrt(gam*gam-1)/gam
      vp0 = gam*vt*sind(alpha)
      vz0 = gam*vt*cosd(alpha)
    case(2) ! Bi-Maxwellian
      read(lun,*) vtperp
      read(lun,*) vtpara
    case(3) 
      !**** Flux map mode ****
      ! This option overrides the 'nparticles' option from
      ! the ptm_parameters.txt file and sets it to the number
      ! of particles necessary to obtain the desired flux map.
      read(lun,*) nEnergies
      read(lun,*) nPitchAngles
      nparticles = nEnergies*nPitchAngles
      write(*,*) 'FLUX MAP MODE'
      write(*,*) 'NPARTICLES = ', nparticles

      allocate(pitchAngles(nPitchAngles))
      allocate(energies(nEnergies))
      read(lun,*) EkevMin
      read(lun,*) EkevMax
      read(lun,*) PitchAngleDegMin
      read(lun,*) PitchAngleDegMax
      if(nEnergies==1) then 
        energies=EkevMin
      else
        energies = linspace(EkevMin,EkevMax,nEnergies)
      endif
      
      if(nPitchAngles==1) then
        pitchAngles = PitchAngleDegMin
      else
        pitchAngles = linspace(PitchAngleDegMin,PitchAngleDegMax,nPitchAngles)
      endif
      
    case default
      call assert(.FALSE.,'particles_initialize','idist='//int2str(idist)//' not supported')
  end select

  close(lun)

  open(unit=lun,file='ptm_input/dist_density_'//id_string//'.txt',action='read',status='old',iostat=ierr)
  call assert(ierr==0,'particles_initialize','error opening dist_density_'//id_string//'.txt')

  read(lun,*) idens
  
  select case(idens)
    case(1) ! All particles start at same position
      read(lun,*) x0
      read(lun,*) y0
      read(lun,*) z0
      write(*,'(a20,3f8.3)') 'x0,y0,z0=', x0, y0, z0
    case(2) ! Particles are seeded randomly throughout a rectangular region
      read(lun,*) xsMin
      read(lun,*) xsMax
      read(lun,*) ysMin
      read(lun,*) ysMax
      read(lun,*) zsMin
      read(lun,*) zsMax
    case(3) ! Particles are randomly seeded at a given radial distance
      read(lun,*) r0
      read(lun,*) mltMin
      read(lun,*) mltMax
    case default
      call assert(.FALSE.,'particles_initialize','idens='//int2str(idens)//' not supported')
  end select

  return

  end subroutine particles_initialize

!

  subroutine particle_cleanup(myParticle)

  implicit none

  type(particle) :: myParticle

  if(associated(myParticle%y)) myParticle%y=>null()

  return

  end subroutine particle_cleanup

!

  subroutine particle_initialize(myParticle,tag)
  ! This is an initialization method for the particle data type
  implicit none

  type(particle) :: myParticle
  integer, intent(in), optional :: tag
  real(dp) :: x, y, z
  real(dp) :: vp, vz, v
  real(dp), dimension(3) :: bvec
  real(dp) :: gam, phi, dx, dy, dz, dt
  integer :: im, jm, km, lm, iEnergy, iPitchAngle

  myParticle%drift = .FALSE.
  myParticle%integrate = .TRUE.
  myParticle%p(1) = e_me*charge_mass_ratio
  myParticle%p(2) = mass*mc2

  select case(idens)
    case(1) ! All particles start at same position
      x = x0
      y = y0
      if(ndim==2) then ! Put the particle in the center of a unit cell in z
        z = 0.5d0
      else ! Full 3d grids are being used
        z = z0
      endif
    case(2) ! Particles are seeded randomly throughout a rectangular region

      x = xsMin+(xsMax-xsMin)*random_uniform()
      y = ysMin+(ysMax-ysMin)*random_uniform()

      if(ndim==2) then
        z = 0.5d0
      else
        z = zsMin+(zsMax-zsMin)*random_uniform()
      endif

    case(3) ! Particles are seeded along the equator at a fixed radial distance and random local times
      phi = 2*pi*(mltMin+random_uniform()*(mltMax-mltMin))/24.0d0
      x = r0*sin(phi)
      y = r0*cos(phi)
      z = merge(0.5d0,0.0d0,ndim==2)
    case default
      call assert(.FALSE.,'particle_initialize','idens='//int2str(idens)//' not supported')
  end select

  select case(idist)
    case(1) ! Initialize a particle beam
      vp = vp0
      vz = vz0
    case(2) ! Bi-Maxwellian
      do ! Get perpendicular momentum is from a truncated normal distribution
        vp = vtperp*random_gauss()
        if(vp>0.d0) exit
      enddo
      ! Parallel momentum comes from a normal Maxwellian
      vz = vtpara*random_gauss()
    case(3) ! Set up particle properties for flux map mode
      call assert(present(tag),'particle_initialize','idist=3 requires specification of particle tag in calling routine')
      iEnergy = ceiling(tag/real(nPitchAngles))
      iPitchAngle = tag-(iEnergy-1)*nPitchAngles
      gam = 1.0+energies(iEnergy)/myParticle%p(2)
      v = ckm*sqrt(gam*gam-1.d0)
      vp = v*sind(pitchAngles(iPitchAngle))
      vz = v*cosd(pitchAngles(iPitchAngle))

      myParticle%fluxMapCoordinates = [merge(thi,tlo,itrace<0),x,y,z,energies(iEnergy),pitchAngles(iPitchAngle)]

    case default
      call assert(.FALSE.,'particle_initialize','idist='//int2str(idist)//' not supported')
  end select

  phi = 2*pi*random_uniform()
  gam = sqrt(1.d0+(vp**2+vz**2)/csq)

  if(ndim==2) vz = 0.d0 ! Just in case of roundoff error

  ! Target component pointers
  allocate(myParticle%y(7))
  myParticle%x => myParticle%y(1:3)
  myParticle%v => myParticle%y(4:6)

  myParticle%x(1) = x
  myParticle%x(2) = y
  myParticle%x(3) = z
  myParticle%v(1) = vp*cos(phi)
  myParticle%v(2) = vp*sin(phi)
  myParticle%v(3) = vz

  ! Set the particle at the correct time in order to allow for forward and reverse tracing
  myParticle%t = merge(Tlo,Thi,itrace>0)
  myParticle%drift = .FALSE.
                       
  ! Initialize the grid and interpolators
  im = locate(xgrid,x)
  jm = locate(ygrid,y)
  km = locate(zgrid,z)
  lm = locate(tgrid,myParticle%t)

  dx = xgrid(im+1)-xgrid(im)
  dy = ygrid(jm+1)-ygrid(jm)
  dz = zgrid(km+1)-zgrid(km)
  dt = tgrid(lm+1)-tgrid(lm)

  call grid_init(myParticle%grid,im,jm,km,lm,dx,dy,dz,dt)

  ! Call the fields routine to initialize the tricubic coefficients and determine the initial timestep
  call get_fields(myParticle,bvec,doInit=.TRUE.)

  if(istep==2) then ! Longer steps with the adaptive RK method
    myParticle%dt = sign(epsilon_drift/abs(myParticle%p(1)*norm(bvec)),real(itrace,dp))
  else
    myParticle%dt = sign(epsilon_orbit/abs(myParticle%p(1)*norm(bvec)),real(itrace,dp))
  endif

  return

  end subroutine particle_initialize

end module particles
