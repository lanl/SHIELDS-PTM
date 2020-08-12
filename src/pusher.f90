!+MODULE NAME
!   stepper
!
!+MODULE DESCRIPTION
!   This module contains routines related to the timestepping of particles. There are
!   routines for both guiding center drifts and full orbit integration and a facility for
!   dynamically switching between these methods.
!
!+AUTHOR
!   Jesse Woodroffe
!   jwoodroffe@lanl.gov
!
!+REVISION HISTORY
!   12 March 2015:    Original code
!    6 April 2015:    Modified code to use RKSUITE90 integrators via local function passing.
!                     As far as I know, only IFORT has this F2008 feature implemented, need to check PGF90
!    5 May 2016:      Modified timestep selection to ensure proper matching of actual times to requested times.
!                     Code has been validated to run on GNU Fortran 4.8.2 with the -std=f2008ts flag
!    3 August 2020:   Heavily streamlined code

module pusher

use global
use rksuite_90
use particles
use fields

contains

  subroutine push(myParticle,tnext)
  ! This is the main particle pushing routine

  implicit none

  type(particle) :: myParticle
  real(dp), intent(in) :: tnext
  real(dp), parameter :: kappa_orbit = 1.0d-2 ! These two values are set so that a particle doesn't just
  real(dp), parameter :: kappa_drift = 0.8d-2 ! oscillate between drift and orbit in regions of marginal smoothness
  real(dp), parameter :: dt_base = 1.0d-6     ! small parameter used to set a maximum number of timesteps
  real(dp) :: kappa, t, t_int
  real(dp), dimension(:), allocatable :: y
  real(dp), dimension(:), allocatable :: yp
  integer :: i, iflag, nstep, nstep_max

  abstract interface
    function deriv_f(t,y) result(ydot)
      integer, parameter :: dp = kind(1.d0)
      real(dp), intent(in) :: t
      real(dp), dimension(:), intent(in) :: y
      real(dp), dimension(size(y)) :: ydot
    end function deriv_f
  end interface

  PROCEDURE(deriv_f), POINTER :: deriv

  ! Initialize storage for the integrators and communicators if using RKSuite (istep=2)

  if(myParticle%drift) then
    if(myParticle%t==merge(Tlo,Thi,itrace>0)) call full_to_drift(myParticle)
    deriv => drift
    allocate(y(4),yp(4))
  else
    deriv => orbit
    allocate(y(6),yp(6))
  endif

  y = myParticle%y(1:size(y))
  if(istep==2) call setup(myParticle%comm,myParticle%t,y,tnext,tol,(/(thresh,i=1,size(y))/),method='M',task='R',message=.FALSE.)

  nstep_max = ceiling(abs(tnext-myParticle%t)/dt_base)

  ! nstep_max is set to be a large fininte number in order to prevent infinite iteration if the timestep approaches zero
  do nstep=1,nstep_max

    t = myParticle%t

	! Part 1 -- Determine Integration Parameters

    kappa = get_kappa(myParticle)

    if (kappa >= kappa_orbit) then ! Use full particle orbits
      if(myParticle%drift) then ! Switch to full orbit
        call drift_to_full(myParticle)
        deriv => orbit
        deallocate(y,yp)
        allocate(y(6),yp(6))
        y = myParticle%y(1:6)
      endif
    else ! Use guiding center equations
      if(.not. myParticle%drift) then ! Switch to guiding center
        call full_to_drift(myParticle)
        deriv => drift
        deallocate(y,yp)
        allocate(y(4),yp(4))
        y = myParticle%y(1:4)
      endif
    endif

    if(istep==2) then ! Modify RKSuite communicator settings
      call collect_garbage(myParticle%comm)
      call setup(myParticle%comm,t,y,tnext,tol,(/(thresh,i=1,size(y))/),method='M',task='R',message=.FALSE.)
    endif

    ! PART 2 -- Timestep Advance

    call calculate_timestep(myParticle)

    ! Check if the timestep is going to take us beyond requested point. If so, set limit of integration appropriately.
    t_int = merge(min(tnext,myParticle%t+myParticle%dt),max(tnext,myParticle%t+myParticle%dt),itrace>0)

    select case(istep)
      case(1) ! Use fixed timestep RK4 integrator
        call rk4(deriv,y,t,t_int-t)
      case(2) ! Use adaptive RKSuite integrator (ADP)
        call range_integrate(myParticle%comm,deriv,t_int,t,y,yp,flag=iflag)
        if(iflag==5) myParticle%integrate = .FALSE.
      case default
        call assert(.FALSE.,"push","istep "//int2str(istep)//" not currently supported.")
    end select

    myParticle%y(1:size(y)) = y
    myParticle%t = t_int

    ! PART 3 -- Boundary Checks

    ! Check if particle is incident on the atmosphere
    if(norm(myParticle%x) <= 1.d0) myParticle%integrate=.FALSE.

    ! Check if particle has reached outer boundaries
    select case(ibound)
      case(1) ! Bounded by limits of domain
        if (myParticle%x(1) < XMin .or. myParticle%x(1) > XMax) myParticle%integrate = .FALSE.
        if (myParticle%x(2) < YMin .or. myParticle%x(2) > YMax) myParticle%integrate = .FALSE.
        if (myParticle%x(3) < ZMin .or. myParticle%x(3) > ZMax) myParticle%integrate = .FALSE.
        if (myParticle%t < TMin .or. myParticle%t > TMax) myParticle%integrate = .FALSE.
      case(2) ! Bounded by plane in negative-x direction
        if (myParticle%x(1) <= xsource) myParticle%integrate=.FALSE.
      case(3) ! Bounded by radial distance
        if (norm(myParticle%x) >= rbound) myParticle%integrate=.FALSE.
      case default
        call assert(.FALSE.,"push","ibound "//int2str(istep)//" not currently supported")
    end select

    if(.not. myParticle%integrate .or. t_int == tnext) exit

  enddo

  ! Clean up integration
  if(istep==2) call collect_garbage(myParticle%comm)
  deallocate(y,yp)

  contains
    ! Making these functions local to the procedure allows us to pass them to the integrator
    ! while maintaining access to the particle information which is necessary to use the
    ! external ODE library RKSuite90. But gfortran does not allow these functions
    ! to be used as arguments in subroutines rk4 and range_integrate above

    function orbit(t,y) result(ydot)
    ! Calculates the derivatives for advancing particle trajectory
    implicit none

    real(dp), dimension(:), intent(in) :: y
    real(dp), dimension(size(y)) :: ydot
    real(dp), intent(in) :: t
    real(dp), dimension(3) :: evec, bvec
    real(dp) :: gami

    if(myParticle%integrate) then

      myParticle%t = t
      myParticle%y(1:6) = y

      call get_fields(myParticle,bvec,evec=evec)
      gami = 1.d0/sqrt(1.d0+dot_product(myParticle%v,myParticle%v)/csq)

      ydot(1:3) =  gami*myParticle%v/re
      ydot(4:6) =  myParticle%qmr*(evec+gami*cross_product(myParticle%v,bvec))

    else
      ! This particle is being taken offline, don't worry about calculating anything.
      ydot = 0.d0

    endif

    return

    end function orbit

    !

    function drift(t,y) result(ydot)
    ! Calculates the derivatives for advancing guiding center trajectory
    ! using the equations of Cary & Brizard [2009]. This version has been modified to inclued
    ! the effects of time-varying magnetic fields via curl E. [method due to M. Henderson]

    implicit none

    real(dp), dimension(:), intent(in) :: y
    real(dp), dimension(size(y)) :: ydot
    real(dp), intent(in) :: t
    real(dp), dimension(3) :: evec, bvec, bhat
    real(dp), dimension(3,3) :: Bgrad, Egrad
    real(dp), dimension(3) :: estar, bstar, gradb, curlb, curlbhat, vdrift, curlE, dbhatdt
    real(dp) :: bstarpara, gami, b0, dupara

    if(myParticle%integrate) then

      myParticle%t = t
      myParticle%y(1:4) = y
      call get_fields(myParticle,bvec,evec=evec,gradb=Bgrad,grade=Egrad)

      b0 = norm(bvec)
      bhat = bvec/b0
      gami = 1.d0/sqrt(1.d0+(2.d0*b0*myParticle%mu+myParticle%upara**2)/csq)

      ! Calculate curls and gradients
      curlb(1) = Bgrad(2,3)-Bgrad(3,2)
      curlb(2) = Bgrad(3,1)-Bgrad(1,3)
      curlb(3) = Bgrad(1,2)-Bgrad(2,1)

      curle(1) = Egrad(2,3)-Egrad(3,2)
      curle(2) = Egrad(3,1)-Egrad(1,3)
      curle(3) = Egrad(1,2)-Egrad(2,1)

      gradb = matmul(Bgrad,bhat)

      curlbhat = (curlb-cross_product(gradb,bhat))/b0
      dbhatdt = (bhat*dot_product(bhat,curle)-curle)/b0

      ! Calculate effective fields
      bstar = bvec+myParticle%upara*curlbhat/myParticle%qmr
      bstarpara = dot_product(bhat,bstar)
      estar = evec-(myParticle%upara*dbhatdt+myParticle%mu*gradb*gami)/myParticle%qmr

      ! Calculate time derivative terms
      vdrift = (gami*myParticle%upara*bstar+cross_product(estar,bhat))/bstarpara
      dupara = myParticle%qmr*dot_product(estar,bstar)/bstarpara

      ydot(1:3) = vdrift/re ! Scale drift velocity to re/s
      ydot(4)   = dupara

    else
      ! This particle has finished integration, take it offline
      ydot = 0.d0

    endif

    return

    end function drift

    !

    subroutine calculate_timestep(myParticle)

    implicit none

    type(particle) :: myParticle
    real(dp) :: usq, gami, b0, wc
    real(dp), dimension(3) :: evec, bvec, bhat, curlb, gradb, curlbhat, fE, fG, fC
    real(dp), dimension(3,3) :: Bgrad

    call get_fields(myParticle,bvec,evec=evec,gradb=Bgrad)
    b0 = norm(bvec)

    if(myParticle%drift) then
      ! Use Kress-Hudson heuristic to pick drift time step
      usq = 2.d0*b0*myParticle%mu+myParticle%upara**2
      gami = 1.d0/sqrt(1.d0+usq/csq)
      bhat = bvec/b0
      curlb(1) = Bgrad(2,3)-Bgrad(3,2)
      curlb(2) = Bgrad(3,1)-Bgrad(1,3)
      curlb(3) = Bgrad(1,2)-Bgrad(2,1)
      gradb = matmul(Bgrad,bhat)
      curlbhat = (curlb-cross_product(gradb,bhat))/b0

      fE = myParticle%qmr*evec                                             ! Electric force
      fG = myParticle%mu*gradb*gami                                        ! Gradient force including mirror term
      fC = myParticle%upara**2*cross_product(curlbhat,bhat)*gami           ! Curvature force

      myParticle%dt = sign(epsilon_drift*sqrt(usq)/norm(fE+fG+fC),real(itrace,dp))
    else
      ! Pick orbit time step to resolve gyro-orbits
      usq = dot_product(myParticle%v,myParticle%v)
      gami = 1.d0/sqrt(1.d0+usq/csq)
      wc = myParticle%qmr*b0*gami
      myParticle%dt = sign(epsilon_orbit/wc,real(itrace,dp))
    endif

    return

    end subroutine calculate_timestep

	!

    function get_kappa(myParticle) result(kappa)

    implicit none

    real(dp) :: kappa
    type(particle), intent(in) :: myParticle
    real(dp) :: b0, usq, uperp, gami, wc, rho
    real(dp), dimension(3) :: bvec, bhat, gradb
    real(dp), dimension(3,3) :: Bgrad

    select case(iswitch)
      case(:-1) ! Guiding center drifts only (GC)
        kappa = 0.9d0*kappa_drift
      case(1:) ! Full particle orbits only (FO)
        kappa = 1.1d0*kappa_orbit
      case(0 ) ! Dynamically switch btwn GC and FO

        call get_fields(myParticle,bvec,gradb=Bgrad)
        b0 = norm(bvec)
        bhat = bvec/b0

        if(myParticle%drift) then
          usq = 2.d0*b0*myParticle%mu+myParticle%upara**2
          uperp = sqrt(2.d0*b0*myParticle%mu)
        else
          usq = dot_product(myParticle%v,myParticle%v)
          uperp = norm(myParticle%v-dot_product(bhat,myParticle%v)*bhat)
        endif

        gami = 1.d0/sqrt(1.d0+usq/csq)
        wc = abs(myParticle%qmr*b0*gami)
        rho = uperp/wc
        gradb = matmul(Bgrad,bhat)
        kappa = rho*norm(gradb)/b0

    end select

    return

    end function get_kappa

  end subroutine push

  !

  subroutine full_to_drift(myParticle)
  ! Switch from full orbit (vx, vy, vz) to guiding center (mu, vpara)

  implicit none

  type(particle) :: myParticle
  real(dp) :: b0, upara, gami, wc
  real(dp), dimension(3) :: bvec, bhat, uperp

  call get_fields(myParticle,bvec)

  ! Determine local magnetic field and it direction
  b0 = norm(bvec)
  bhat = bvec/b0
  gami = 1.d0/sqrt(1.d0+dot_product(myParticle%v,myParticle%v)/csq)

  wc = gami*myParticle%qmr*b0

  upara= dot_product(bhat,myParticle%v)
  uperp = myParticle%v-upara*bhat

  ! Move the particle to the guiding center
  myParticle%x = myParticle%x + cross_product(myParticle%v,bhat)/(re*wc)

  myParticle%upara => myParticle%v(1)
  myParticle%mu => myParticle%v(2)

  ! Get magnetic field at new location
  call get_fields(myParticle,bvec)
  b0 = norm(bvec)

  ! Determine current drift parameters
  myParticle%upara = upara
  myParticle%mu = dot_product(uperp,uperp)/(2.d0*b0)

  myParticle%drift = .TRUE.

  return

  end subroutine full_to_drift

!

  subroutine drift_to_full(myParticle)
  ! This version corrects an error that arose because the sign function takes the absolute value
  ! of its first argument; this caused the velocity and position vectors to be not orthogonal and
  ! may have caused errors in trajectories.

  implicit none

  type(particle) :: myParticle
  real(dp), dimension(3) :: rhat, vhat, bvec, bhat, xgc, rkeep
  real(dp), dimension(3) :: gradb
  real(dp), dimension(3,3) :: R, Bgrad
  real(dp) :: vperp, vpara, rho, wc, b0, phi, dbmin, dbnow, gami
  integer :: i

  call get_fields(myParticle,bvec)
  b0 = norm(bvec)
  bhat = bvec/b0
  gami = 1.d0/sqrt(1.d0+2.d0*b0*myParticle%mu+myParticle%upara**2)

  wc = abs(myParticle%qmr*b0*gami)

  vpara = myParticle%upara
  vperp = sqrt(2.d0*b0*myParticle%mu)
  rho = vperp/(re*wc)

  select case(iphase)

    case(1) ! Random gyrophase
      phi = twopi*random_uniform()

      ! Get gyrovector in Cartesian coordinates
      R = fac_rotation_matrix(myParticle,toCartesian=.TRUE.)
      rhat = matmul(R,[cos(phi),sin(phi),0.0d0])

      vhat = sign(1.d0,myParticle%qmr)*cross_product(bhat,rhat)
      vhat = vhat/norm(vhat)

    case(2) ! Brute-force minimization of delta B

      R = fac_rotation_matrix(myParticle,toCartesian=.TRUE.)
      xgc = myParticle%x

      do i=1,nphase
        phi = twopi*real(i,dp)
        rhat = [cos(phi),sin(phi),0.d0]
        rhat = matmul(R,rhat)
        myParticle%x = xgc+rho*rhat

        call get_fields(myParticle,bvec)
        dbnow = abs(b0-norm(bvec))

        if(i==1) then
          dbmin = dbnow
          rkeep = rhat
        else
          if(dbnow < dbmin) then
            dbmin = dbnow
            rkeep = rhat
          endif
        endif

      enddo

      rhat = rkeep
      vhat = sign(1.d0,myParticle%qmr)*cross_product(rhat,bhat)
      vhat = vhat/norm(vhat)
      myParticle%x = xgc

  case(3) ! Pfefferle's grad B method

    ! Calculate the local magnetic field and the direction of its gradient
    call get_fields(myParticle,bvec,gradb=Bgrad)

    gradb = matmul(Bgrad,bhat)

    ! Put particle along direction orthogonal to both gradb and b
    rhat = cross_product(bhat,gradb)
    rhat = rhat/norm(rhat)

    vhat = (gradb-bhat*dot_product(bhat,gradb))/norm(cross_product(bvec,gradb))
    vhat = sign(1.d0,myParticle%qmr)*cross_product(rhat,bhat)
    vhat = vhat/norm(vhat)

  end select

  myParticle%x = myParticle%x+rho*rhat
  myParticle%v = vpara*bhat+vperp*vhat
  myParticle%drift = .FALSE.

  if(associated(myParticle%upara)) then
    myParticle%upara => NULL()
  endif

  if(associated(myParticle%mu)) then
    myParticle%mu => NULL()
  endif

  return

  end subroutine drift_to_full

  !

  subroutine rk4(derivs,y,t,dt)
  ! Basic 4th order Runge-Kutta for timestepping. This is fixed-step with no error estimation.
  ! An adaptive method with error control is available from the RKSuite library (set istep=2)
  !
  ! Jesse Woodroffe, 8/7/15

  implicit none

  interface
    function derivs(t,y)
    use global, only: dp
    real(dp), intent(in) :: t
    real(dp), dimension(:), intent(in) :: y
    real(dp), dimension(size(y)) :: derivs
    end function derivs
  end interface

  real(dp), dimension(:), intent(inout) :: y
  real(dp), intent(in) :: t
  real(dp), intent(in) :: dt
  real(dp), dimension(size(y)) :: k1, k2, k3, k4

  k1 = dt*derivs(t,y)
  k2 = dt*derivs(t+0.5d0*dt,y+0.5d0*k1)
  k3 = dt*derivs(t+0.5d0*dt,y+0.5d0*k2)
  k4 = dt*derivs(t+dt,y+k3)

  y = y+(k1+2.d0*k2+2.d0*k3+k4)/6.d0

  return

  end subroutine rk4

end module pusher
