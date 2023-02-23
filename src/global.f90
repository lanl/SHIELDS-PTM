!+MODULE NAME
!   global
!
!+MODULE DESCRIPTION
!   This module contains global constants and variables that need to be
!   accessed by multiple other modules along with miscellaneous utility
!   routines.
!
!+AUTHOR
!   Jesse Woodroffe
!   jwoodroffe@lanl.gov
!
!+REVISION HISTORY
!   12    March 2015:   Original code
!   31    March 2015:   Revised code to use new particle type
!    4      May 2015:   Revised particle object for convenient linear interpolation in time
!    4 February 2016:   Replaced signal_error error handler with assert

module global

use rksuite_90

integer, parameter :: dp = kind(1.d0)

real(dp), parameter :: ckm = 2.998d5   ! speed of light in km/s
real(dp), parameter :: csq = ckm*ckm
real(dp), parameter :: mc2 = 5.11d2    ! electron rest mass energy in keV
real(dp), parameter :: pi = 4.d0*atan(1.d0)
real(dp), parameter :: dtor = pi/180.d0    ! degrees to radians
real(dp), parameter :: e_me = 1.758d2    ! electron charge to mass ratio mult. by 1e9
real(dp), parameter :: me_e = 1.d0/e_me    ! electron mass to charge ratio (x 1e9)
real(dp), parameter :: re = 6371.d0    ! Volumetric Earth radius in km
real(dp), parameter :: tol = 1.0d-6
real(dp), parameter :: thresh = 1.0d-10
real(dp), parameter :: twopi = 8.d0*atan(1.d0)
real(dp), parameter :: epsilon_drift = 0.1d0    ! factor for setting timestepping in GC approximation
real(dp), parameter :: epsilon_orbit = 0.05d0    ! factor for setting timestepping in full orbit
real(dp), dimension(3,3), parameter :: eye = reshape([1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0],[3,3])

integer :: runid
integer :: iseed
integer :: ndim
integer :: istep
integer :: iswitch
integer :: iphase
integer :: nphase
integer :: itrace
integer :: iFirst
integer :: iLast
integer :: ntot
integer :: nx
integer :: ny
integer :: nz
integer :: nt
integer :: nparticles
integer :: itraj
integer :: ibound

character(len=4) :: id_string     ! why here ? move to PTM.f90 subroutine get_run_id

real(dp) :: dtIn
real(dp) :: dtOut
real(dp) :: XMin
real(dp) :: XMax
real(dp) :: YMin
real(dp) :: YMax
real(dp) :: ZMin
real(dp) :: ZMax
real(dp) :: TMin
real(dp) :: TMax
real(dp) :: charge
real(dp) :: mass
real(dp) :: charge_mass_ratio
real(dp) :: tlo
real(dp) :: thi
real(dp) :: xsource
real(dp) :: rbound

logical :: fluxMap = .FALSE.

! Gyrophase calculation parameters
real(dp) :: dphi

! Spatial and temporal grids
real(dp), dimension(:), allocatable :: xgrid
real(dp), dimension(:), allocatable :: ygrid
real(dp), dimension(:), allocatable :: zgrid
real(dp), dimension(:), allocatable :: tgrid

! Four-dimensional (3 space, 1 time) EM Fields
real(dp), dimension(:,:,:,:), allocatable :: BX3D
real(dp), dimension(:,:,:,:), allocatable :: BY3D
real(dp), dimension(:,:,:,:), allocatable :: BZ3D
real(dp), dimension(:,:,:,:), allocatable :: EX3D
real(dp), dimension(:,:,:,:), allocatable :: EY3D
real(dp), dimension(:,:,:,:), allocatable :: EZ3D

! OBJECT TYPE DEFINITIONS

type grid
  ! This object stores all of the information about the local grid cell. The grid is a
  ! sub-object of the particle
  integer :: im   ! Lower x-index
  integer :: jm   ! Lower y-index
  integer :: km   ! lower z-index
  integer :: lm   ! lower t-index
  real(dp) :: dx  ! Size of grid cell in x
  real(dp) :: dy  ! Size of grid cell in y
  real(dp) :: dz  ! Size of grid cell in z
  real(dp) :: dt  ! Size of grid cell in t
end type grid

type particle
  real(dp) :: t                                             ! Current time for particle
  real(dp) :: dt                                            ! Timestep for integrator
  integer :: tag                                            ! Particle number
  logical :: drift,integrate                                ! Logical switches
  real(dp), dimension(:), pointer :: x                      ! Position
  real(dp), dimension(:), pointer :: v                      ! Relativistic velocity
  real(dp), dimension(:), pointer :: y                      ! Actual storage
  real(dp), pointer :: mu, upara, qmr, mc2
  real(dp), dimension(:), pointer:: p                       ! Parameters: p(1) = charge/mass ratio; p(2) = particle mass in keV
  real(dp), dimension(2,64) :: bxinterp, byinterp, bzinterp ! Tricubic interpolation coefficients for B
  real(dp), dimension(2,64) :: exinterp, eyinterp, ezinterp ! Tricubic interpolation coefficients for E
  type(rk_comm_real_1d) :: comm                             ! Communicator for use with RKSuite
  type(grid) :: grid                                        ! 4D grid cell used by interpolators
  real(dp), dimension(9) :: fluxMapCoordinates              ! Liouville coordinates for quick flux map
end type particle

interface read_array
  ! These are routines to read 1-, 2-, and 3-dimensional binary data files
  ! with data stored in C-style.
  module procedure read_array_1d, read_array_2d, read_array_3d
end interface read_array

contains

  subroutine assert(condition,caller,message)
  ! Logical assertion-based error handling. To signal a failure
  ! without any test, just call assert with condition=.FALSE.
  implicit none

  logical, intent(in) :: condition
  character(len=*) :: caller, message
  integer :: nrep

  if(.not. condition) then

    nrep = max(len(caller)+22,len(message)+19)

    write(*,*)
    write(*,*) repeat('_',nrep)
    write(*,*) repeat('-',nrep)
    write(*,*) "Error encountered in: "//caller
    write(*,*) "Error message was: "//message
    write(*,*) "Program terminated by ASSERT"
    write(*,*) repeat('_',nrep)
    write(*,*) repeat('-',nrep)
    write(*,*)

    error stop 1

  endif

  return

  end subroutine assert

!

  function int2str(i) result(istr)
  ! Create a string version of an integer number
  implicit none

  integer, intent(in) :: i
  character(len=1+floor(log10(real(i,dp)))) :: istr
  character(len=8) :: ilong

  call assert(1+floor(log10(real(i,dp))) <= 8,'int2str','long integer data type not supported')

  write(ilong,'(i8)') i

  istr = trim(adjustl(ilong))

  return

  end function int2str

!++++++++++++++++++++++++++++++++++++++++
! INTERFACE READ_ARRAY PROCEDURES
!++++++++++++++++++++++++++++++++++++++++

  subroutine read_array_3d(fname,dat3d)
  ! Data is stored in C-order, so we have to reorder the data in order to get
  ! it in a shape we understand. Note that this uses "STREAM" access, a feature
  ! of the Fortran 2003 standard.
  implicit none

  character(len=*), intent(in) :: fname
  real(dp), dimension(:,:,:), intent(out) :: dat3d
  real(dp), dimension(size(dat3d,3),size(dat3d,2),size(dat3d,1)) :: temp
  integer :: j, lun

  open(newunit=lun,file=fname,access='stream',form='unformatted')
  read(lun) temp
  close(lun)

  do j=1,size(dat3d,2)
      dat3d(:,j,:) = transpose(temp(:,j,:))
  enddo

  return

  end subroutine read_array_3d

  !

  subroutine read_array_2d(fname,dat2d)

  implicit none

  character(len=*), intent(in) :: fname
  real(dp), dimension(:,:), intent(out) :: dat2d
  real(dp), dimension(size(dat2d,1),size(dat2d,2)) :: temp
  integer :: lun

  open(newunit=lun,file=fname,access='stream',form='unformatted')
  read(lun) temp
  close(lun)

  dat2d = temp

  return

  end subroutine read_array_2d

  !

  subroutine read_array_1d(fname,dat1d)

  implicit none

  character(len=*), intent(in) :: fname
  real(dp), dimension(:), intent(out) :: dat1d
  integer :: lun

  open(newunit=lun,file=fname,access='stream',form='unformatted')
  read(lun) dat1d
  close(lun)

  return

  end subroutine read_array_1d

!++++++++++++++++++++++++++++++++++++++++
! Reading procedure for arrays having
! unspecified length
!++++++++++++++++++++++++++++++++++++++++

  subroutine read_grid(fname,grid)
  ! This routine determines the number of lines in an ascii file,
  ! allocates an array for the appropriate number of values, then
  ! reads in the values into the newly-allocated array.

  ! Use the iso_fortran_env module to get the end-of-file status value
  use iso_fortran_env

  implicit none

  character(len=*), intent(in) :: fname
  real(dp), dimension(:), allocatable, intent(out) :: grid
  integer :: lun, i, istat, n

  open(file=fname,newunit=lun,status='old',iostat=istat)
  call assert(istat==0,'read_grid','Error opening '//fname)

  ! Count the number of lines in the file
  n=0
  do
    read(unit=lun,fmt=*,iostat=istat)
    if(istat == iostat_end) exit
    call assert(istat==0,'read_grid','Error reading from '//fname)
    n=n+1
  enddo

  ! Allocate the array with appropriate number of values
  allocate(grid(n))

  ! Return to the beginning of the file and read it into the array
  rewind(lun)
  read(lun,*) grid
  close(lun)

  return

  end subroutine read_grid

!++++++++++++++++++++++++++++++++++++++++
! Utility Functions
!++++++++++++++++++++++++++++++++++++++++

  function norm(v)
  ! Determine the magnitude of a vector
  implicit none

  real(dp), dimension(3) :: v
  real(dp) :: norm

  norm = sqrt(dot_product(v,v))

  return

  end function norm

!

  function cross_product(a,b) result(c)
  ! Calculate the cross product of two vectors
  implicit none

  real(dp), dimension(3), intent(in) :: a,b
  real(dp), dimension(3) :: c

  c(1) = a(2)*b(3)-a(3)*b(2)
  c(2) = a(3)*b(1)-a(1)*b(3)
  c(3) = a(1)*b(2)-a(2)*b(1)

  return

  end function cross_product

!

  function outer_product(a,b) result(M)
  ! Calculate the outer product of two vectors
  implicit none

  real(dp), dimension(:), intent(in) :: a,b
  real(dp), dimension(size(a),size(b)) :: M
  integer :: i, j

  do i=1,size(a)
    do j=1,size(b)
      M(i,j) = a(i)*b(j)
    enddo
  enddo

  return

  end function outer_product

!

  function linspace(a,b,n) result(v)
  ! Find a linearly spaced vector of length n between a and b.
  ! This function copies the functionality of Matlab's linspace routine.

  real(dp), intent(in) :: a, b
  integer, intent(in) :: n
  real(dp), dimension(n) :: v
  real(dp) :: dv
  integer :: i

  call assert(n>0,'linspace','length of requested vector = '//int2str(n))

  if(n>1) then
    dv = (b-a)/real(n-1,dp)
    do i=1,n
      v(i) = a+real(i-1,dp)*dv
    enddo
  else
    call assert(a/=b,'linspace','different endpoints for single element vector')
    v(1) = a
  endif

  return

  end function linspace

!++++++++++++++++++++++++++++++++++++++++
! Angular conversion functions
!++++++++++++++++++++++++++++++++++++++++

  function degrees(rads) result(v)

  real(dp), intent(in) :: rads
  real(dp) :: v

  v = rads/dtor

  end function degrees

!

  function radians(degs) result(v)

  real(dp), intent(in) :: degs
  real(dp) :: v

  v = degs*dtor

  end function radians

!++++++++++++++++++++++++++++++++++++++++
! Degree-argument trigonometric functions
!++++++++++++++++++++++++++++++++++++++++

  function sind(x)

  real(dp), intent(in) :: x
  real(dp) :: sind

  sind = sin(dtor*x)

  return

  end function sind

!

  function cosd(x)

  implicit none

  real(dp), intent(in) :: x
  real(dp) :: cosd

  cosd = cos(dtor*x)

  end function cosd

!++++++++++++++++++++++++++++++++++++++++
! Random Number Generators
!++++++++++++++++++++++++++++++++++++++++

  function random_uniform() result(dev)
  ! This is a functionalized wrapper of Fortran's intrinsic
  ! random number generator.
  implicit none

  real(dp) :: dev

  call random_number(dev)

  return

  end function random_uniform

!

  function random_gauss() result(harvest)
  ! Gaussian random deviate using Box-Mueller Method (NRF90)
  implicit none

  real(dp) :: harvest
  real(dp) :: rsq,v1,v2
  real(dp), SAVE :: g
  logical :: stored=.false.

  if (stored) then
    harvest=g
    stored=.false.
  else
    do
      call random_number(v1)
      call random_number(v2)
      v1=2.0d0*v1-1.0d0
      v2=2.0d0*v2-1.0d0
      rsq=v1**2+v2**2
      if (rsq > 0.0d0 .and. rsq < 1.0d0) exit
    end do
    rsq=sqrt(-2.0d0*log(rsq)/rsq)
    harvest=v1*rsq
    g=v2*rsq
    stored=.true.
  end if

  end function random_gauss

end module global
