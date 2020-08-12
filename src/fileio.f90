!+MODULE NAME
!   fileio
!
!+MODULE DESCRIPTION
!   This module contains routines associated with reading and writing data files
!
!+AUTHOR
!   Jesse Woodroffe
!   jwoodroffe@lanl.gov

module fileio

use global
use fields
use particles

contains

  subroutine read_ptm_parameters()

  implicit none

  character(len=4) :: test_string
  integer :: lun, ierr

  open(newunit=lun,file='ptm_input/ptm_parameters_'//id_string//'.txt',action='read',status='old',iostat=ierr)

  call assert(ierr==0,'read_ptm_parameters','could not open ptm_parameters_'//id_string//'.txt')

  read(lun,*) runid
  read(lun,*) iseed
  read(lun,*) nparticles
  read(lun,*) ndim
  read(lun,*) itrace
  read(lun,*) ifirst
  read(lun,*) ilast
  read(lun,*) ntot
  read(lun,*) dtIn
  read(lun,*) dtOut
  read(lun,*) istep
  read(lun,*) iswitch
  read(lun,*) iphase
  read(lun,*) nphase
  read(lun,*) charge
  read(lun,*) mass
  read(lun,*) tlo
  read(lun,*) thi
  read(lun,*) itraj
  read(lun,*) ibound
  if(ibound == 2) then
    read(lun,*) xsource
  else if(ibound==3) then
    read(lun,*) rbound
  endif

  ! Gyrophase switching. dphi is defined in the GLOBAL module
  call assert(nphase > 0,'read_ptm_parameters','nphase must be greater than zero.')
  dphi = twopi/real(nphase,dp)

  ! Charge/mass ratio for all particles
  charge_mass_ratio = charge/mass

  ! Number of timesteps
  nt = ilast-ifirst+1

  close(lun)

  write(*,*) 'RUNID = ', runid

  write(test_string,'(i4.4)') runid

  call assert(test_string==id_string,'read_ptm_parameters','runid inconsistent between filename and file contents')

  return

  end subroutine read_ptm_parameters

!

  subroutine storeData(myParticle,dataStore)
  ! Determine output data and put it in a structure for later writing
  implicit none

  type(particle) :: myParticle
  real(dp), dimension(8), intent(out) :: dataStore
  real(dp), dimension(3) :: bvec, bhat, xpos, gradb, rvec
  real(dp), dimension(3,3) :: Bgrad
  real(dp) :: b0, upara, uperp, gam

  call get_fields(myParticle,bvec,gradb=Bgrad)
  b0 = norm(bvec)
  bhat = bvec/b0

! Store lab frame quantities = orbit params

  if(myParticle%t==merge(Tlo,Thi,itrace>0)) then
    upara = dot_product(bhat,myParticle%v)
    uperp = norm(myParticle%v-upara*bhat)
    gam=sqrt(1.d0+(upara/ckm)**2+(uperp/ckm)**2)
  else
    if(myParticle%drift) then
      upara = myParticle%upara
      uperp = sqrt(2.0d0*b0*myParticle%mu)
    else
      upara = dot_product(bhat,myParticle%v)
      uperp = norm(myParticle%v-upara*bhat)
    endif
  endif

  gam=sqrt(1.d0+(upara/ckm)**2+(uperp/ckm)**2)

  dataStore(1) = myParticle%t
  dataStore(2:4) = myParticle%x
  dataStore(5) = uperp/gam
  dataStore(6) = upara/gam
  dataStore(7) = myParticle%mc2*(gam-1.d0)
  dataStore(8) = degrees(atan2(uperp,upara))

  return

  end subroutine storeData

!

  subroutine writeDataStore(dataStore,tag)
  ! This routine is called by a single thread at a time, so the SAVED firstCall variable is okay: the first thread
  ! to call it should set its status for all other subsequent calls. Theoretically, two threads might try to access
  ! this routine simultaneously (and that could cause problems), so we have to make calls to this routine inside
  ! an OMP CRITICAL structure.

  implicit none

  real(dp), dimension(:,:), intent(in) :: dataStore
  integer, intent(in) :: tag
  integer :: i, lun, ierr
  logical, save :: firstCall = .TRUE.

  if(firstCall) then
    open(newunit=lun,file='ptm_output/ptm_'//id_string//'.dat',action='write',status='replace')
    firstCall = .FALSE.
  else
    open(newunit=lun,file='ptm_output/ptm_'//id_string//'.dat',action='write',status='old',position='append',iostat=ierr)
    call assert(ierr==0,'writeDataStore','error opening ptm_'//int2str(runid)//'.dat')
  endif

  write(lun,*) '# ', tag
  ! The # symbol separates data from different particles

  do i=lbound(dataStore,1),ubound(dataStore,1)
    write(lun,'(8es17.7e3)') dataStore(i,:)
  enddo

  close(lun)

  return

  end subroutine writeDataStore

!

  subroutine writeFluxCoordinates(myParticle)
  ! This program writes out a simple set of flux mapping coordinates that can be used to apply Liouville's theorem

  implicit none

  type(particle), intent(in) :: myParticle
  integer :: lun, ierr
  logical, save :: firstCall = .TRUE.
  real(dp) :: gam, b0, energy
  real(dp), dimension(3) :: bvec

  if(firstCall) then
    open(newunit=lun,file='ptm_output/map_'//id_string//'.dat',action='write',status='replace')
    write(lun,'(4es17.7e3)') myParticle%fluxMapCoordinates(1:4)
    firstCall = .FALSE.
  else
     open(newunit=lun,file='ptm_output/map_'//id_string//'.dat',action='write',status='old',position='append',iostat=ierr)
     call assert(ierr==0,'writeFluxCoordinates','error opening ptm_'//int2str(runid)//'.dat')
  endif

  if(myParticle%drift) then
    call get_fields(myParticle,bvec)
    b0 = norm(bvec)
    gam=sqrt(1+(myParticle%upara/ckm)**2+2*b0*myParticle%mu/csq)
  else
    gam=sqrt(1+dot_product(myParticle%v,myParticle%v)/csq)
  endif

  energy = (gam-1)*myParticle%mc2

  ! Write fluxmap data. Time (1), final position (2,3,4), ?? (5,6),
  ! ?final/initial? energy (7), initial velocity vector (8,9,10)
  write(lun,'(10es17.7e3)') myParticle%t, myParticle%x, myParticle%fluxMapCoordinates(5:6), energy, myParticle%fluxMapCoordinates(7:9)

  close(lun)

  end subroutine writeFluxCoordinates

end module fileio
