module finite_differences
!+Description
!  This module contains routines to approximate mixed first derivatives to a discretely sampled function
!  using repeated one-dimensional parabolic interpolation. The function need not be sampled uniformly, but
!  it is assumed that the grid is Cartesian with coordinate grids that are independent of all other coordinates
!  (e.g. this does not work for AMR/nested grids)
!
!  The purpose of these routines is to approximate derivatives for use with the tricubic interpolation routines
!  that are used in our particle stepper. This hard-coded 2nd order method is much faster than the N-th order
!  Lagrange method. This is very similar to the finite_differences_qp module that had been previously used,
!  but a somewhat different formulation of the coefficients is used and it is a little more flexible.
!
!+Routines
!  deriv1    Approximate the first derivative of a 1-d function 
!  deriv2    Approximate the mixed 2nd derivative of a 2-d function
!  deriv3    Approximate the mixed 3rd derivative of a 3-d function
!
!+Author
!  Jesse Woodroffe
!  jwoodroffe@lanl.gov

use global

contains

  function deriv1(x,f,idx) result(der)
  !+Description 
  !  Approximates the derivatives of a discretely sampled function by fitting a parabola 
  !  and then differentiating the parabola:
  !
  !   f =  A(x-x0)^2 + B(x-x0) + C
  !  fp = 2A(x-x0) + B
  !
  !  The optional parameter idx allows us to specify if the derivative at the left (idx < 0) or right (idx>0) 
  !  cell is desired instead of the center value. Derivatives at the center are slightly more accurate, but this
  !  option allows us to get derivatives at the edge of the grid. 
  !
  !+Inputs
  !  x(3) double precision array of grid points
  !  f(3) double precision array of function values at the grid points
  !  idx  integer that specifies the grid point at which the derivative is to be evaluated
  !
  !+Outputs
  !  der  double precision scalar with derivative at desired grid point
  !
  ! Author:
  !  Jesse Woodroffe
  !  jwoodroffe@lanl.gov

  implicit none

  real(dp), dimension(3) :: x
  real(dp), dimension(3) :: f
  real(dp) :: A, B, p, q, dxp, dxm, dxtot, dfp, dfm, denom
  real(dp) :: der
  integer, intent(in) :: idx

  dxp = x(3)-x(2)
  dxm = x(2)-x(1)
  dxtot = dxp+dxm
  denom = 1.d0/(dxtot*dxp*dxm)

  dfp = f(3)-f(2)
  dfm = f(1)-f(2)

  p = dfm*dxp
  q = dfp*dxm

  B = denom*(dxm*q-dxp*p)

  der = 0.d0
  if(idx==0) then ! Evaluate at center
    der = B
  else if(idx==1) then ! Evaluate at right
    A = denom*(p+q)
    der = B+2*A*dxp
  else if(idx==-1) then ! Evaluate at left
    A = denom*(p+q)
    der = B-2*A*dxm
  endif
    
  return

  end function deriv1

!

  function deriv2(x,y,f,idx,idy) result(der)
  !+Description
  !  Calculate 2-d derivative using repeated application of 1d derivative.
  !
  !+Inputs
  !  x(3)   double precision array of grid points for dimension 1
  !  y(3)   double precision array of grid points for dimension 2
  !  f(3,3) double precision array of function values at the grid points
  !  idx    integer that specifies the grid point at which the derivative is to be evaluated in x
  !  idy    integer that specifies the grid point at which the derivative is to be evaluated in y
  !
  !+Outputs
  !  der  double precision scalar with derivative at desired grid point specified by idx, idy
  !
  !+Author
  !  Jesse Woodroffe
  !  jwoodroffe@lanl.gov
  !
  implicit none

  real(dp), dimension(3), intent(in) :: x, y
  real(dp), dimension(3,3), intent(in) :: f
  real(dp), dimension(3) :: temp  
  real(dp) :: der
  integer, intent(in) :: idx, idy
  integer :: i 
  
  do i=1,3
    temp(i) = deriv1(x,f(:,i),idx=idx)
  enddo

  der = deriv1(y,temp,idx=idy)

  return

  end function deriv2

!

  function deriv3(x,y,z,f,idx,idy,idz) result(der)
  !+Description
  !  Calculate 3-d derivative using repeated application of 1d derivative.
  !
  !+Inputs
  !  x(3)     double precision array of grid points for dimension 1
  !  y(3)     double precision array of grid points for dimension 2
  !  z(3)     double precision array of grid points for dimension 3
  !  f(3,3,3) double precision array of function values at the grid points
  !  idx      integer that specifies the grid point at which the derivative is to be evaluated in x
  !  idy      integer that specifies the grid point at which the derivative is to be evaluated in y
  !  idz      integer that specifies the grid point at which the derivative is to be evaluated in z
  !
  !+Outputs
  !  der  double precision scalar with derivative at desired grid point specified by idx, idy, idz
  !
  !+Author
  !  Jesse Woodroffe
  !  jwoodroffe@lanl.gov
  !
  implicit none

  real(dp), intent(in), dimension(3) :: x, y, z
  real(dp), intent(in), dimension(3,3,3) :: f
  real(dp), dimension(3) :: temp
  real(dp) :: der
  integer, intent(in) :: idx, idy, idz

  integer :: i

  do i=1,3
    temp(i) = deriv2(y,z,f(i,:,:),idx=idy,idy=idz)
  enddo

  der = deriv1(x,temp,idx=idx)

  return

  end function deriv3

end module finite_differences
