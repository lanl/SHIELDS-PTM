!+MODULE NAME
!   particles
!
!+MODULE DESCRIPTION
!   This module initialize the particle spatial and velocity distributions,
!   and contains a function for transforming particle location between Cartesian
!   and field-aligned coordinate systems.
!
!+CONTAINS subroutines
!   PARTICLES_INITIALIZE
!   PARTICLE_INITIALIZE
!   FAC_ROTATION_MATRIX
!
!+AUTHOR
!   Jesse Woodroffe
!   jwoodroffe@lanl.gov
!

module particles

use global
use rksuite_90
use fields
use interpolation

integer, private :: idist, idens
real(dp), private :: x0, y0, z0, E0, alpha, vtperp, vtpara, vz0, vp0
real(dp), private :: xsMin, xsMax, ysMin, ysMax, zsMin, zsMax, mltMin, mltMax, r0

! Flux mapping
integer, private :: nEnergies, nPitchAngles
real(dp), private :: EkevMin, EkevMax, PitchAngleDegMin, PitchAngleDegMax
real(dp), dimension(:), allocatable, private :: pitchAngles, energies
real(dp), private :: phi0 ! Initial gyrophase angle in degrees; if phi0 < 0, use random gyrophase.

contains

  subroutine particles_initialize()

  implicit none

  integer :: nseed, lun, ierr
  real(dp) :: gam, vt
  integer :: i

  ! Configure velocity space distribution

  open(newunit=lun,file='ptm_input/dist_velocity_'//id_string//'.txt',action='read',status='old',iostat=ierr)
  call assert(ierr==0,'particles_initialize','error opening dist_velocity_'//id_string//'.txt')

  ! Initialize the system Random Number Generator
  call random_seed(size=nseed)
  call random_seed(put=(/(iseed+i,i=0,nseed-1)/))

  ! Read in the user-specified distribution function type
  read(lun,*) idist

  select case(idist)

    case(1)
      ! Monoenergetic ring distribution
      read(lun,*) E0
      read(lun,*) alpha
      read(lun,*) phi0
      gam = 1.d0+E0/(mass*mc2)
      vt = ckm*sqrt(gam*gam-1)/gam
      vp0 = vt*sind(alpha)
      vz0 = vt*cosd(alpha)

    case(2)
      ! Non-relativistic Bi-Maxwellian distribution
      read(lun,*) vtperp
      read(lun,*) vtpara
      read(lun,*) phi0

    case(3)
      ! Automatic flux map mode
      ! This option overrides the 'nparticles' option from
      ! the ptm_parameters.txt file and sets it to the number
      ! of particles necessary to obtain the desired flux map.

      fluxMap = .TRUE.

      read(lun,*) nEnergies
      read(lun,*) nPitchAngles
      read(lun,*) phi0
      read(lun,*) EkevMin
      read(lun,*) EkevMax
      read(lun,*) PitchAngleDegMin
      read(lun,*) PitchAngleDegMax
      read(lun,*) xsource

      nparticles = nEnergies*nPitchAngles
      write(*,*) 'FLUX MAP MODE'
      write(*,*) 'NPARTICLES = ', nparticles

      allocate(pitchAngles(abs(nPitchAngles)))
      allocate(energies(abs(nEnergies)))

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

    case(4)
      ! User-specified flux map mode
      ! This option overrides the 'nparticles' option from
      ! the ptm_parameters.txt file and sets it to the number
      ! of particles necessary to obtain the desired flux map.

      fluxMap = .TRUE.
      read(lun,*) nEnergies
      read(lun,*) nPitchAngles
      read(lun,*) phi0
      read(lun,*) xsource

      nparticles = nEnergies*nPitchAngles
      write(*,*) 'FLUX MAP MODE'
      write(*,*) 'NPARTICLES = ', nparticles

      allocate(pitchAngles(abs(nPitchAngles)))
      allocate(energies(abs(nEnergies)))
      call read_array('ptm_data/energies.bin',energies)
      call read_array('ptm_data/pitchangles.bin',pitchAngles)

    case default

      call assert(.FALSE.,'particles_initialize','idist='//int2str(idist)//' not supported')

  end select

  close(lun)


  ! Configure physical space distributions

  open(newunit=lun,file='ptm_input/dist_density_'//id_string//'.txt',action='read',status='old',iostat=ierr)
  call assert(ierr==0,'particles_initialize','error opening dist_density_'//id_string//'.txt')

  read(lun,*) idens

  select case(idens)
    case(1) ! All particles are seeded at same position
      read(lun,*) x0
      read(lun,*) y0
      read(lun,*) z0
      write(*,'(a20,3f8.3)') 'x0,y0,z0=', x0, y0, z0
    case(2) ! Particles are seeded randomly throughout a cubic region
      read(lun,*) xsMin
      read(lun,*) xsMax
      read(lun,*) ysMin
      read(lun,*) ysMax
      read(lun,*) zsMin
      read(lun,*) zsMax
    case(3) ! Particles are seeded randomly at a given radial distance
      read(lun,*) r0
      read(lun,*) mltMin
      read(lun,*) mltMax
    case default
      call assert(.FALSE.,'particles_initialize','idens='//int2str(idens)//' not supported')

  end select

  close(lun)

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

  function fac_rotation_matrix(myParticle,toCartesian) result(R)
  ! Calculate a matrix that rotates from Cartesian ("XYZ") coordinates to a between a local
  ! field-aligned coordinate system ("PQB"). PQB is equaivalent to the dipole coordinate system
  ! when the magnetic field is dipolar. The inverse transform is given by the transpose of the
  ! matrix.
  !
  ! P is the "x-like" coordinate which points in a direction perpendicular to both Q and B
  !   directions. For a dipole field, this is in the direction of the local radius of curvature
  ! Q is the "y-like" coordinate which points in a direction perpendicular to both the radial
  !   direction and the magnetic field. For a dipole magnetic field this is the azimuthal direction.
  ! B is the "z-like" coordinate which points in the local magnetic field direction
  !
  ! Jesse Woodroffe
  ! 8/22/2016

  implicit none

  type(particle) :: myParticle
  real(dp), dimension(3,3) :: R
  real(dp), dimension(3) :: bhat, phat, qhat, rhat
  logical, intent(in), optional :: toCartesian

  call get_fields(myParticle,bhat)
  bhat = bhat/norm(bhat)
  rhat = myParticle%x/norm(myParticle%x)
  qhat = cross_product(bhat,rhat)
  phat = cross_product(qhat,bhat)

  ! Populate rotation matrix
  R(1,:) = phat
  R(2,:) = qhat
  R(3,:) = bhat

  ! If rotation from FAC to Cartesian is desired, transpose rotation matrix
  if(present(toCartesian)) then
    if(toCartesian) R = transpose(R)
  endif

  return

  end function fac_rotation_matrix

!

  subroutine particle_initialize(myParticle,tag)
  ! This is an initialization method for the particle data type
  implicit none

  type(particle) :: myParticle
  integer, intent(in), optional :: tag
  real(dp) :: x, y, z
  real(dp) :: vp, vz, v
  real(dp), dimension(3) :: bvec, bhat, rhat
  real(dp), dimension(3,3) :: R
  real(dp) :: gam, phi, dx, dy, dz, dt, wc
  integer :: im, jm, km, lm, iEnergy, iPitchAngle

  real(dp) :: usq, gami, b0
  real(dp), dimension(3) :: evec, curlb, curlbhat, gradb, fE, fG, fC
  real(dp), dimension(3,3) :: Bgrad

  myParticle%drift = .FALSE.
  myParticle%integrate = .TRUE.

  allocate(myParticle%p(2))
  allocate(myParticle%y(7))
  myParticle%x => myParticle%y(1:3)
  myParticle%v => myParticle%y(4:6)

  myParticle%p(1) = e_me*charge_mass_ratio
  myParticle%p(2) = mass*mc2

  myParticle%qmr => myParticle%p(1)
  myParticle%mc2 => myParticle%p(2)

! initialize particle location (x,y,z)

  select case(idens)

    case(1) ! All particles start at same position
      x = x0
      y = y0
      z = z0
    case(2) ! Particles are seeded randomly throughout a rectangular region
      x = xsMin+(xsMax-xsMin)*random_uniform()
      y = ysMin+(ysMax-ysMin)*random_uniform()
      z = zsMin+(zsMax-zsMin)*random_uniform()
    case(3) ! Particles are seeded along the equator at a fixed radial distance and random local times
      phi = 2*pi*(mltMin+random_uniform()*(mltMax-mltMin))/24.0d0
      x = r0*sin(phi)
      y = r0*cos(phi)
      z = 0.0d0
    case default
      call assert(.FALSE.,'particle_initialize','idens='//int2str(idens)//' not supported')
  end select

  myParticle%x = [x,y,z]

  ! initialize particle velocity components (vp,vz)

  select case(idist)
    case(1) ! Initialize a particle beam
      vp = vp0
      vz = vz0
      gam = 1.0/sqrt(1.d0-(vp**2+vz**2)/csq)

    case(2) ! Bi-Maxwellian
      do ! Get perpendicular velocity truncated normal distribution
        vp = vtperp*random_gauss()
        if(vp>0.d0) exit
      enddo
      ! Parallel velocity from a normal distribution
      vz = vtpara*random_gauss()
      gam = 1.0/sqrt(1.d0-(vp**2+vz**2)/csq)

    case(3:4) ! Set up particle properties for flux map mode
      call assert(present(tag),'particle_initialize','idist=3 requires specification of particle tag in calling routine')
      myParticle%tag = tag
      iEnergy = ceiling(tag/real(nPitchAngles))
      iPitchAngle = tag-(iEnergy-1)*nPitchAngles
      gam = 1.0+energies(iEnergy)/myParticle%mc2
      v = ckm*sqrt(gam*gam-1.d0)/gam
      vp = v*sind(pitchAngles(iPitchAngle))
      vz = v*cosd(pitchAngles(iPitchAngle))
      myParticle%fluxMapCoordinates(1:6) = [merge(thi,tlo,itrace<0),x,y,z,energies(iEnergy),pitchAngles(iPitchAngle)]

    case default
      call assert(.FALSE.,'particle_initialize','idist='//int2str(idist)//' not supported')

  end select

  ! *********** construct myParticle object

  ! Set the particle at the correct time in order to allow for forward and reverse tracing
  myParticle%t = merge(Tlo,Thi,itrace>0)
  myParticle%drift = iswitch<0

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
  call get_fields(myParticle,bvec,doInit=.TRUE.)

  bhat = bvec/norm(bvec)                                  ! mag field direction
  R = fac_rotation_matrix(myParticle,toCartesian=.TRUE.)  ! rotation from PQB field-aligned to Cartesian coordinates
  phi = merge(radians(phi0),twopi*random_uniform(),phi0>0.d0)
  rhat = matmul(R,[cos(phi),sin(phi),0.d0])

  myParticle%v = gam*(vz*bhat+vp*rhat)

  ! Store initial cartesian components for flux mapping
  myParticle%fluxMapCoordinates(7:9) = myParticle%v

  return

  end subroutine particle_initialize

end module particles
