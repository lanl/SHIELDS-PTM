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
use stepper
use omp_lib

implicit none

integer :: n, iwrite, nwrite
real(dp) :: tic, toc
type(particle) :: myParticle
real(dp), dimension(:,:), allocatable :: particleData

call get_run_id()
call read_ptm_parameters()
call fields_initialize() 
call particles_initialize()
call initialize_timing()

! Zero index for initial conditions
allocate(particleData(0:nwrite,8))

!$omp parallel do default(shared) private(n,myParticle,particleData,iwrite)
do n=1,nparticles

  call particle_initialize(myParticle,n)
  call storeData(myParticle,particleData(0,:)) ! Initial conditions
  
  do iwrite=1,nwrite
    call stepper_push(myParticle,myParticle%t+sign(dtOut,real(itrace,dp)))
    call storeData(myParticle,particleData(iwrite,:))
    if(.not. myParticle%integrate) exit
  enddo
!$omp critical
  if(iwrite < nwrite) then 
    ! Particle encountered a boundary and kicked out early
    call writeDataStore(particleData(:iwrite,:))
  else 
    ! Successful execution of all timesteps
    call writeDataStore(particleData)
  endif
  call writeFluxCoordinates(myParticle)
!$omp end critical
  call particle_cleanup(myParticle)
enddo
!$omp end parallel do

deallocate(particleData)

call fields_cleanup()
call finalize_timing()

write(*,*) "PTM simulation has finished"

contains

  subroutine get_run_id()

  implicit none

  integer :: argnum, arglen, argstat
  character(len=3) :: argstr
  argnum = 1
  arglen = 3
  call get_command_argument(argnum,argstr,arglen,argstat)

  id_string = repeat('0',4-len(trim(argstr)))//trim(argstr)

  return

  end subroutine get_run_id

  !

  subroutine initialize_timing()

  implicit none

  nwrite = ceiling((THi-TLo)/abs(dtOut))

  write(*,*) "Fields Time Range = ", TMin, TMax
  write(*,*) "Particle Time Range = ", TLo, THi
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
