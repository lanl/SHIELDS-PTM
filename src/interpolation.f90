!+MODULE NAME
!  INTERPOLATION
!
!+MODULE DESCRIPTION
!   This module contains routines for trilinear and tricubic interpolation of 3D gridded
!   data. Trilinear interpolation is a triply-iterated linear interpolation. Tricubic is
!   an implementation of a local tricubic interpolation scheme based on the paper
!   by F. Leiken and J. Marsden, "Tricubic Interpolation in Three Dimensions".
!   This implementation differs from others (e.g., iterated cubic spline) because it is purely local -- 
!   rather than using multiple points, it relies on multiple derivatives. When these derivatives are
!   not explicitly available, they can be approximated using finite difference methods. Currently,
!   finite differences are approximated using the FINITE_DIFFERENCES_CS2 module which estimates the 
!   derivatives using a quadratic polynomial fit.
!
!   This code is not currently configured to directly handle AMR grids, so such data must be
!   pre-processed onto a cartesian mesh.
!
!+GLOBAL DATA
!  GRID data type
!
!+ROUTINES
!  TRICUBIC_INIT          Subroutine    -> Calcluate tricubic interpolation coefficients for a given grid cell
!  GRID_INIT              Subroutine    -> Initialize a grid object
!  TRILINEAR_INTERPOLATE  Subroutine    -> Evaluate a trilinear interpolant using the data in a GRID object
!  TRICUBIC_INTERPOLATE   Subroutine    -> Evaluate a tricubic interpolant using pre-calculated coefficients from TRICUBIC_INIT and the data in a GRID object
!  LOCATE                 Function      -> Search an ordered array for the lower bracketing index, i.e. the index of the
!                                          element of the array whose value is closest to a specified value without going over
!
!+DEPENDENCIES
!  NONE
!
!+USAGE
!   1. Initialize the grid by calling GRID_INIT()
!   2. Determine the local interpolation coefficients by calling TRICUBIC_INIT() 
!      with appropriate inputs
!   3. Calculate the local value of the desired quantities using trilinear_interpolate() and
!      tricubic_interpolate() routines with appropriate inputs
!   4. Note that repeated evaluations of the tricubic interpolant can be carried out for different values 
!      of the independent variables so long as they remain within the region of validity for the grid object
!
!+AUTHOR
!  Jesse Woodroffe
!  jwoodroffe@lanl.gov
!
!+REVISION HISTORY
!   16 January 2015: Original code
!   17 January 2015: Added tricubic_init and removed calculation of B_INV from tricubic_setup
!   18 January 2015: Added matrix_inverse module with matrix inversion routine 
!                    Added INTERPOLANT data type and code snippets for use in future implementations
!   20 January 2015: Changed code to fully support use of interpolant data type and added test program
!   21 January 2015: ::CODE FORKED TO VERSION 2:: Revised derivative scaling in test program. Modified 
!                    interpolant object and moved derivative scaling to the evaluation procedures so 
!                    that interpolation is transparent to the user. The example program has been 
!                    modified to minimize the number of times tricubic_setup() has to be called.
!   22 January 2015: ::CODE FORKED TO VERSION 3:: Added setup and evaluation methods to INTERPOLANT class.
!                    The interpolation process has been converted to a more object-oriented flow that
!                    permits the simultaneous existence of multiple interpolants. Added an explicit definiton
!                    of the B_INV matrix and removed the routines originally used to calculate it.
!   23 January 2015: Added locate function to search ordered but non-uniformly spaced ("stretched") grids.       
!   11 March 2015:   Fixed an error in some of the mixed derivative routines
!   17 March 2015:   Intoroduced GRID object to streamline data sharing amongst interpolants. Added trilinear
!                    interpolation routine for electric fields
!   04 March 2015:   Changed interpolant object to 64-element array, altered some variable names.
!   13 August 2015:  Fixed error in B_INV array

module interpolation

use global

! The B_INV array is a bunch of constants that are used to determine the interpolation weights.
! There is no need to recalculate this array, so I've included it here as a monstrously large constant.
! Note that the square bracket [] notation for array definition is only supported by the FORTRAN standard
! since 2008, but the Intel ifort compiler doesn't complain.

! The B_INV matrix relates the derivatives to the coefficients of the interpolant. This is a geometrical constant for a 1 x 1 x 1 cube, so no need to actually calculate it on the fly.
! The bracket assignment notation is not supported by all Fortran compilers, so use with others may require us to put the array in a file and read it once at the beginning of the simulation.
real(dp),dimension(64,64), parameter :: B_INV = reshape([&
[   1,   0,  -3,   2,   0,   0,   0,   0,  -3,   0,   9,  -6,   2,   0,  -6,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   0,   9,  -6,   0,   0,   0,   0,   9,   0, -27,  18,  -6,   0,  18, -12,   2,   0,  -6,   4,   0,   0,   0,   0,  -6,   0,  18, -12,   4,   0, -12,   8], &
[   0,   0,   3,  -2,   0,   0,   0,   0,   0,   0,  -9,   6,   0,   0,   6,  -4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -9,   6,   0,   0,   0,   0,   0,   0,  27, -18,   0,   0, -18,  12,   0,   0,   6,  -4,   0,   0,   0,   0,   0,   0, -18,  12,   0,   0,  12,  -8], &
[   0,   0,   0,   0,   0,   0,   0,   0,   3,   0,  -9,   6,  -2,   0,   6,  -4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -9,   0,  27, -18,   6,   0, -18,  12,   0,   0,   0,   0,   0,   0,   0,   0,   6,   0, -18,  12,  -4,   0,  12,  -8], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   9,  -6,   0,   0,  -6,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, -27,  18,   0,   0,  18, -12,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  18, -12,   0,   0, -12,   8], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,   0,  -9,   6,   0,   0,   0,   0,  -9,   0,  27, -18,   6,   0, -18,  12,  -2,   0,   6,  -4,   0,   0,   0,   0,   6,   0, -18,  12,  -4,   0,  12,  -8], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   9,  -6,   0,   0,   0,   0,   0,   0, -27,  18,   0,   0,  18, -12,   0,   0,  -6,   4,   0,   0,   0,   0,   0,   0,  18, -12,   0,   0, -12,   8], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   9,   0, -27,  18,  -6,   0,  18, -12,   0,   0,   0,   0,   0,   0,   0,   0,  -6,   0,  18, -12,   4,   0, -12,   8], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  27, -18,   0,   0, -18,  12,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, -18,  12,   0,   0,  12,  -8], &
[   0,   1,  -2,   1,   0,   0,   0,   0,   0,  -3,   6,  -3,   0,   2,  -4,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   6,  -3,   0,   0,   0,   0,   0,   9, -18,   9,   0,  -6,  12,  -6,   0,   2,  -4,   2,   0,   0,   0,   0,   0,  -6,  12,  -6,   0,   4,  -8,   4], &
[   0,   0,  -1,   1,   0,   0,   0,   0,   0,   0,   3,  -3,   0,   0,  -2,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,  -3,   0,   0,   0,   0,   0,   0,  -9,   9,   0,   0,   6,  -6,   0,   0,  -2,   2,   0,   0,   0,   0,   0,   0,   6,  -6,   0,   0,  -4,   4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,  -6,   3,   0,  -2,   4,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -9,  18,  -9,   0,   6, -12,   6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   6, -12,   6,   0,  -4,   8,  -4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   3,   0,   0,   2,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   9,  -9,   0,   0,  -6,   6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -6,   6,   0,   0,   4,  -4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,  -6,   3,   0,   0,   0,   0,   0,  -9,  18,  -9,   0,   6, -12,   6,   0,  -2,   4,  -2,   0,   0,   0,   0,   0,   6, -12,   6,   0,  -4,   8,  -4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   3,   0,   0,   0,   0,   0,   0,   9,  -9,   0,   0,  -6,   6,   0,   0,   2,  -2,   0,   0,   0,   0,   0,   0,  -6,   6,   0,   0,   4,  -4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   9, -18,   9,   0,  -6,  12,  -6,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -6,  12,  -6,   0,   4,  -8,   4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -9,   9,   0,   0,   6,  -6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   6,  -6,   0,   0,  -4,   4], &
[   0,   0,   0,   0,   1,   0,  -3,   2,  -2,   0,   6,  -4,   1,   0,  -3,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   0,   9,  -6,   6,   0, -18,  12,  -3,   0,   9,  -6,   0,   0,   0,   0,   2,   0,  -6,   4,  -4,   0,  12,  -8,   2,   0,  -6,   4], &
[   0,   0,   0,   0,   0,   0,   3,  -2,   0,   0,  -6,   4,   0,   0,   3,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -9,   6,   0,   0,  18, -12,   0,   0,  -9,   6,   0,   0,   0,   0,   0,   0,   6,  -4,   0,   0, -12,   8,   0,   0,   6,  -4], &
[   0,   0,   0,   0,   0,   0,   0,   0,  -1,   0,   3,  -2,   1,   0,  -3,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,   0,  -9,   6,  -3,   0,   9,  -6,   0,   0,   0,   0,   0,   0,   0,   0,  -2,   0,   6,  -4,   2,   0,  -6,   4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   2,   0,   0,   3,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   9,  -6,   0,   0,  -9,   6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -6,   4,   0,   0,   6,  -4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,   0,  -9,   6,  -6,   0,  18, -12,   3,   0,  -9,   6,   0,   0,   0,   0,  -2,   0,   6,  -4,   4,   0, -12,   8,  -2,   0,   6,  -4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   9,  -6,   0,   0, -18,  12,   0,   0,   9,  -6,   0,   0,   0,   0,   0,   0,  -6,   4,   0,   0,  12,  -8,   0,   0,  -6,   4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   0,   9,  -6,   3,   0,  -9,   6,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,  -6,   4,  -2,   0,   6,  -4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -9,   6,   0,   0,   9,  -6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   6,  -4,   0,   0,  -6,   4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,  -3,   2,   0,   0,   0,   0,  -3,   0,   9,  -6,   2,   0,  -6,   4,  -2,   0,   6,  -4,   0,   0,   0,   0,   6,   0, -18,  12,  -4,   0,  12,  -8,   1,   0,  -3,   2,   0,   0,   0,   0,  -3,   0,   9,  -6,   2,   0,  -6,   4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,  -2,   0,   0,   0,   0,   0,   0,  -9,   6,   0,   0,   6,  -4,   0,   0,  -6,   4,   0,   0,   0,   0,   0,   0,  18, -12,   0,   0, -12,   8,   0,   0,   3,  -2,   0,   0,   0,   0,   0,   0,  -9,   6,   0,   0,   6,  -4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,   0,  -9,   6,  -2,   0,   6,  -4,   0,   0,   0,   0,   0,   0,   0,   0,  -6,   0,  18, -12,   4,   0, -12,   8,   0,   0,   0,   0,   0,   0,   0,   0,   3,   0,  -9,   6,  -2,   0,   6,  -4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   9,  -6,   0,   0,  -6,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, -18,  12,   0,   0,  12,  -8,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   9,  -6,   0,   0,  -6,   4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   0,   3,  -2,   0,   0,   0,   0,   3,   0,  -9,   6,  -2,   0,   6,  -4,   1,   0,  -3,   2,   0,   0,   0,   0,  -3,   0,   9,  -6,   2,   0,  -6,   4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   2,   0,   0,   0,   0,   0,   0,   9,  -6,   0,   0,  -6,   4,   0,   0,   3,  -2,   0,   0,   0,   0,   0,   0,  -9,   6,   0,   0,   6,  -4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   0,   9,  -6,   2,   0,  -6,   4,   0,   0,   0,   0,   0,   0,   0,   0,   3,   0,  -9,   6,  -2,   0,   6,  -4], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -9,   6,   0,   0,   6,  -4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   9,  -6,   0,   0,  -6,   4], &
[   0,   0,   0,   0,   0,   1,  -2,   1,   0,  -2,   4,  -2,   0,   1,  -2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   6,  -3,   0,   6, -12,   6,   0,  -3,   6,  -3,   0,   0,   0,   0,   0,   2,  -4,   2,   0,  -4,   8,  -4,   0,   2,  -4,   2], &
[   0,   0,   0,   0,   0,   0,  -1,   1,   0,   0,   2,  -2,   0,   0,  -1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,  -3,   0,   0,  -6,   6,   0,   0,   3,  -3,   0,   0,   0,   0,   0,   0,  -2,   2,   0,   0,   4,  -4,   0,   0,  -2,   2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   2,  -1,   0,   1,  -2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,  -6,   3,   0,  -3,   6,  -3,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -2,   4,  -2,   0,   2,  -4,   2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,  -1,   0,   0,  -1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   3,   0,   0,   3,  -3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,  -2,   0,   0,  -2,   2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,  -6,   3,   0,  -6,  12,  -6,   0,   3,  -6,   3,   0,   0,   0,   0,   0,  -2,   4,  -2,   0,   4,  -8,   4,   0,  -2,   4,  -2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   3,   0,   0,   6,  -6,   0,   0,  -3,   3,   0,   0,   0,   0,   0,   0,   2,  -2,   0,   0,  -4,   4,   0,   0,   2,  -2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   6,  -3,   0,   3,  -6,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,  -4,   2,   0,  -2,   4,  -2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,  -3,   0,   0,  -3,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -2,   2,   0,   0,   2,  -2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,  -2,   1,   0,   0,   0,   0,   0,  -3,   6,  -3,   0,   2,  -4,   2,   0,  -2,   4,  -2,   0,   0,   0,   0,   0,   6, -12,   6,   0,  -4,   8,  -4,   0,   1,  -2,   1,   0,   0,   0,   0,   0,  -3,   6,  -3,   0,   2,  -4,   2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   1,   0,   0,   0,   0,   0,   0,   3,  -3,   0,   0,  -2,   2,   0,   0,   2,  -2,   0,   0,   0,   0,   0,   0,  -6,   6,   0,   0,   4,  -4,   0,   0,  -1,   1,   0,   0,   0,   0,   0,   0,   3,  -3,   0,   0,  -2,   2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,  -6,   3,   0,  -2,   4,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -6,  12,  -6,   0,   4,  -8,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,  -6,   3,   0,  -2,   4,  -2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   3,   0,   0,   2,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   6,  -6,   0,   0,  -4,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   3,   0,   0,   2,  -2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   2,  -1,   0,   0,   0,   0,   0,   3,  -6,   3,   0,  -2,   4,  -2,   0,   1,  -2,   1,   0,   0,   0,   0,   0,  -3,   6,  -3,   0,   2,  -4,   2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,  -1,   0,   0,   0,   0,   0,   0,  -3,   3,   0,   0,   2,  -2,   0,   0,  -1,   1,   0,   0,   0,   0,   0,   0,   3,  -3,   0,   0,  -2,   2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   6,  -3,   0,   2,  -4,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,  -6,   3,   0,  -2,   4,  -2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,  -3,   0,   0,  -2,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   3,   0,   0,   2,  -2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,  -3,   2,  -2,   0,   6,  -4,   1,   0,  -3,   2,   0,   0,   0,   0,  -2,   0,   6,  -4,   4,   0, -12,   8,  -2,   0,   6,  -4,   0,   0,   0,   0,   1,   0,  -3,   2,  -2,   0,   6,  -4,   1,   0,  -3,   2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,  -2,   0,   0,  -6,   4,   0,   0,   3,  -2,   0,   0,   0,   0,   0,   0,  -6,   4,   0,   0,  12,  -8,   0,   0,  -6,   4,   0,   0,   0,   0,   0,   0,   3,  -2,   0,   0,  -6,   4,   0,   0,   3,  -2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   0,   3,  -2,   1,   0,  -3,   2,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,  -6,   4,  -2,   0,   6,  -4,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   0,   3,  -2,   1,   0,  -3,   2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   2,   0,   0,   3,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   6,  -4,   0,   0,  -6,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   2,   0,   0,   3,  -2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   0,   3,  -2,   2,   0,  -6,   4,  -1,   0,   3,  -2,   0,   0,   0,   0,   1,   0,  -3,   2,  -2,   0,   6,  -4,   1,   0,  -3,   2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   2,   0,   0,   6,  -4,   0,   0,  -3,   2,   0,   0,   0,   0,   0,   0,   3,  -2,   0,   0,  -6,   4,   0,   0,   3,  -2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,  -3,   2,  -1,   0,   3,  -2,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   0,   3,  -2,   1,   0,  -3,   2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   3,  -2,   0,   0,  -3,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   2,   0,   0,   3,  -2], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,  -2,   1,   0,  -2,   4,  -2,   0,   1,  -2,   1,   0,   0,   0,   0,   0,  -2,   4,  -2,   0,   4,  -8,   4,   0,  -2,   4,  -2,   0,   0,   0,   0,   0,   1,  -2,   1,   0,  -2,   4,  -2,   0,   1,  -2,   1], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   1,   0,   0,   2,  -2,   0,   0,  -1,   1,   0,   0,   0,   0,   0,   0,   2,  -2,   0,   0,  -4,   4,   0,   0,   2,  -2,   0,   0,   0,   0,   0,   0,  -1,   1,   0,   0,   2,  -2,   0,   0,  -1,   1], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   2,  -1,   0,   1,  -2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,  -4,   2,   0,  -2,   4,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   2,  -1,   0,   1,  -2,   1], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,  -1,   0,   0,  -1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -2,   2,   0,   0,   2,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,  -1,   0,   0,  -1,   1], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   2,  -1,   0,   2,  -4,   2,   0,  -1,   2,  -1,   0,   0,   0,   0,   0,   1,  -2,   1,   0,  -2,   4,  -2,   0,   1,  -2,   1], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,  -1,   0,   0,  -2,   2,   0,   0,   1,  -1,   0,   0,   0,   0,   0,   0,  -1,   1,   0,   0,   2,  -2,   0,   0,  -1,   1], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,  -2,   1,   0,  -1,   2,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   2,  -1,   0,   1,  -2,   1], &
[   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   1,   0,   0,   1,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,  -1,   0,   0,  -1,   1]],[64,64])

contains

  subroutine grid_init(myGrid, im, jm, km, lm, dx, dy, dz, dt)
  ! Initialize a grid object used in trilinear and tricubic interpolation
  implicit none

  type(grid) :: myGrid
  integer, intent(in) :: im, jm, km, lm
  real(dp), intent(in) :: dx, dy, dz, dt

  myGrid%im = im
  myGrid%jm = jm
  myGrid%km = km
  myGrid%lm = lm
  myGrid%dx = dx
  myGrid%dy = dy
  myGrid%dz = dz
  myGrid%dt = dt

  return

  end subroutine grid_init

!

  function trilinear_interpolate(myGrid,xw,yw,zw,F) result(h)
  ! Find the value of a three-dimensional function at a point using
  ! trilinear interpolation of gridded values
  implicit none

  type(grid) :: myGrid
  real(dp), intent(in) :: xw, yw, zw
  real(dp), dimension(2,2,2) :: F
  real(dp) :: x, y, z
  real(dp) :: f00, f10, f01, f11
  real(dp) :: g0, g1
  real(dp) :: h

  x = (xw-xgrid(myGrid%im))/myGrid%dx
  y = (yw-ygrid(myGrid%jm))/myGrid%dy
  z = (zw-zgrid(myGrid%km))/myGrid%dz

  f00 = (1-x)*F(1,1,1)+x*F(2,1,1)
  f10 = (1-x)*F(1,2,1)+x*F(2,2,1)
  f01 = (1-x)*F(1,1,2)+x*F(2,1,2)
  f11 = (1-x)*F(1,2,2)+x*F(2,2,2)
  g0  = (1-y)*f00+y*f10
  g1  = (1-y)*f01+y*f11
  h   = (1-z)*g0+z*g1

  return

  end function trilinear_interpolate

!

  subroutine tricubic_init(myCoeffs,myGrid,F,FX,FY,FZ,FXY,FXZ,FYZ,FXYZ)
  ! Calculate the interpolation coefficients for an interpolant object
  implicit none
  
  real(dp), dimension(64) :: myCoeffs
  type(grid) :: myGrid
  real(dp), dimension(2,2,2), intent(in) :: F,FX,FY,FZ,FXY,FXZ,FYZ,FXYZ
  real(dp), dimension(64) :: bbar
  integer :: ix, iy, iz, irow
       
  !Determine the coefficients of the interpolation using abar = B_INV . bbar
  ! We only need to save the results of this calculation, bbar is not needed later.

  irow=0
  do iz=1,2
    do iy=1,2
      do ix=1,2
        irow=irow+1
        bbar(irow   ) =    F(ix,iy,iz)
        bbar(irow+8 ) =   FX(ix,iy,iz)*myGrid%dx
        bbar(irow+16) =   FY(ix,iy,iz)*myGrid%dy
        bbar(irow+24) =   FZ(ix,iy,iz)*myGrid%dz
        bbar(irow+32) =  FXY(ix,iy,iz)*myGrid%dx*myGrid%dy
        bbar(irow+40) =  FXZ(ix,iy,iz)*myGrid%dx*myGrid%dz
        bbar(irow+48) =  FYZ(ix,iy,iz)*myGrid%dy*myGrid%dz
        bbar(irow+56) = FXYZ(ix,iy,iz)*myGrid%dx*myGrid%dy*myGrid%dz
      enddo
    enddo
  enddo

  ! This is a matrix-vector multiply between a 64x64 matrix and a 64 element vector. Tests indicate that
  ! vendor-optimized matmuls (such as mkl) are more efficient than a call to the L2 BLAS DGEMV routine
  ! (results may differ for an optimized BLAS library, but the times weren't even particularly close).

  myCoeffs = matmul(B_INV,bbar)
  
  return
  
  end subroutine tricubic_init

!

   subroutine tricubic_interpolate(myCoeffs, myGrid, x0, y0, z0, f0, gradf)
  ! Evaluate a tricubic interpolant at a point. You can also request the gradient of
  ! the function by including the optional gradf vector.
  implicit none

  real(dp), dimension(64) :: myCoeffs
  type(grid) :: myGrid
  real(dp), intent(in) :: x0, y0, z0
  real(dp), intent(out) :: f0
  real(dp), dimension(3), intent(out), optional :: gradf
  real(dp) :: x, y, z
  real(dp) :: xpow, ypow, zpow
  real(dp) :: xder, yder, zder
  integer :: i, j, k, idex
  
  ! Map the simulation cell to a unit cube
  x = (x0-xgrid(myGrid%im))/myGrid%dx
  y = (y0-ygrid(myGrid%jm))/myGrid%dy
  z = (z0-zgrid(myGrid%km))/myGrid%dz

  ! Calculate the value of the function based on tricubic interpolation coefficients
  ! previously calculated during a call to tricubic_setup. Also calculate the gradient
  ! of the function if the gradf argument is present

  f0 = 0.d0
  if(present(gradf)) gradf = 0.d0

  ! These #pow and #der variables are used to avoid using the exponentiation operator, 
  ! and they turn out to be a fairly efficient way of expressing the algorithm.

  idex=0
  zder = 0.d0
  zpow = 1.d0
  do k=0,3
    yder = 0.d0
    ypow = 1.d0
    do j=0,3
      xder = 0.d0
      xpow = 1.d0
      do i=0,3
        idex = idex+1
        f0 = f0+myCoeffs(idex)*xpow*ypow*zpow
        if(present(gradf)) then ! Calculate the gradient
          gradf(1) = gradf(1)+myCoeffs(idex)*i*xder*ypow*zpow
          gradf(2) = gradf(2)+myCoeffs(idex)*j*xpow*yder*zpow
          gradf(3) = gradf(3)+myCoeffs(idex)*k*xpow*ypow*zder
        endif
        xder = xpow
        xpow = xpow*x
      enddo
      yder = ypow
      ypow = ypow*y
    enddo
    zder = zpow
    zpow = zpow*z
  enddo
  
  if(present(gradf)) then ! Scale derivatives by actual size of the cell
    gradf(1) = gradf(1)/myGrid%dx
    gradf(2) = gradf(2)/myGrid%dy
    gradf(3) = gradf(3)/myGrid%dz
  endif

  return
  
  end subroutine tricubic_interpolate

!

  function locate(v,x,ih) result(i0)
  ! Seach for the lower bracketing index, i.e. the index of array element of a monotonically
  ! increasing array v that is closest to the value x without being larger than x. If you know
  ! approximately where to look, you can specify an optional "hunting index". This function
  ! is essentially an amalgam of the Numerical Recipes "locate" and "hunt" routines.
  implicit none

  real(dp), dimension(:), intent(in) :: v
  real(dp), intent(in) :: x
  integer, intent(in), optional :: ih
  integer :: dh
  integer :: i0
  integer :: ilo, ihi, imid
  integer :: n  

  n = size(v)

  ! Check that the value is actually bracketed. It's up to the user to handle these exceptional cases. 
  if(x<v(1)) then
    i0 = 0
    return
  else if(x>v(n)) then
    i0 = n+1      
    return
  endif

  ! If the user has specified a starting index, we will hunt from that point to determine the bracket
  if(present(ih)) then ! Hunt
    dh = 1 ! Initialize the search range to one cell; if we're lucky, we'll already be bracketing the value
    ihi = ih
    do
      ilo = ihi
      ihi = ilo+dh
      if(ihi > n) then ! That's as far as we can go, truncate search range to last element
        ihi = n
        exit
      endif
      if(v(ihi) > x) exit ! Hunt is successful
      dh = 2*dh ! The hunt hasn't yet succeeded, so we expand the search range and try again
    enddo
    ilo = ih
    ihi = ih+dh
  else ! We don't know where to start, so use the entire interval
    ilo = 1
    ihi = n
  endif
  ! Bisection search
  do
    if(ihi - ilo == 1) exit
    imid = (ihi+ilo)/2
    if(v(imid) == x) then
      ilo = imid
      exit
    else if(v(imid) > x) then
      ihi = imid
    else
      ilo = imid
    endif
  enddo
  i0 = ilo
  return
  end function locate
  
end module interpolation
