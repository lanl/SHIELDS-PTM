!+PROGRAM NAME
!   particle_test
!
!+PROGRAM DESCRIPTION
!   particle_test is a driver routine for the SHIELDS-PTM code
!
!+AUTHOR
!   Jesse Woodroffe
!   jwoodroffe@lanl.gov
!

program particle_test

use global
use fileio
use particles
use fields
use pusher
use omp_lib

implicit none

integer :: n, iwrite, nwrite
real(dp) :: tic, toc, lastepoch
type(particle) :: myParticle
real(dp), dimension(:,:), allocatable :: particleData

call get_run_id()

call read_ptm_parameters()

call fields_initialize()            ! reads (t,x,y,z)-grid and initial/final B and E fields from ptm_data files

call particles_initialize()         ! reads bounds from ptm_input files, needed for particle init below (do loop)

call initialize_timing()

allocate(particleData(0:nwrite,8))  ! zero index for initial conditions, eight output qunatities of myParticle: t,x,y,z,vperp,vpara,E,PA


!$omp parallel do default(shared) private(n,myParticle,particleData,iwrite)

DO n=1,nparticles

  call particle_initialize(myParticle,n)  ! initializes velocity (vapara,vperp) and spatial distribution in object myParticle

  call storeData(myParticle,particleData(0,:)) ! store initial conditions

  do iwrite=1,nwrite
    if (mod(iwrite,360).eq.0) write (*,*) 'particle ',n,' time =',myParticle%t
    call push(myParticle,myParticle%t+sign(dtOut,real(itrace,dp)))   !  advance myParticle one dtOut (in file ptm_pars)
    call storeData(myParticle,particleData(iwrite,:))                        !  store data every output file cadence
    ! exit when finished or when problems were encountered with particle time-stepping in stepper_push
    if(.not. myParticle%integrate .or. iwrite==nwrite) exit
    lastepoch=myParticle%t
  enddo
!$omp critical

  ! write data files
  if(fluxMap) call writeFluxCoordinates(myParticle)
  if((.not. fluxMap) .or. itraj==1) then
    call writeDataStore(particleData(0:iwrite,:),n)
  end if

!$omp end critical

  call particle_cleanup(myParticle)

END DO

!$omp end parallel do

write(*,*) "PTM simulation has finished at time t =",lastepoch
write(*,*) "Desired end time =", max(Tlo,itrace*Thi)

deallocate(particleData)
call fields_cleanup()
call finalize_timing()


contains

  subroutine get_run_id()

  implicit none

  integer :: argnum, arglen, argstat
  character(len=4) :: argstr            ! set at run
! character(len=4) :: id_string         ! output of this subr, not declared here, but in GLOBAL.f90

  argnum = 1
  arglen = 4
  call get_command_argument(argnum,argstr,arglen,argstat)

  id_string = repeat('0',4-len(trim(argstr)))//trim(argstr)

  return

  end subroutine get_run_id

  !

  subroutine initialize_timing()

  implicit none

  nwrite = ceiling((THi-TLo)/abs(dtOut))+1

  write(*,*) "Fields Set at Time Range = ", TMin, TMax
  write(*,*) "Simulation Time Range = ", TLo, THi
  write(*,*) "dtOut = ", dtOut
  write(*,*) "nwrite = ", nwrite
  write(*,*) "Done with setup, beginning PTM simulation"

  tic = omp_get_wtime()

  end subroutine initialize_timing

  !

  subroutine finalize_timing()

  implicit none

  toc = omp_get_wtime()

  write(*,*) "PTM simulation duration = ", toc-tic, " seconds"

  return

  end subroutine finalize_timing

end program particle_test
