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
!    6 April 2015:    Modified code to use RKSUITE90 integrators via local function passing. As far as I know,
!                     only IFORT has this F2008 feature implemented, need to check PGF90
!    5 May 2016:      Modified timestep selection to ensure proper matching of actual times to requested times. Code
!                     has been validated to run on GNU Fortran 4.8.2 with the -std=f2008ts flag

module stepper

use global
use rksuite_90
use particles
use fields

contains

  subroutine stepper_push(myParticle,tnext)
  ! This is the main particle pushing routine. Calculates the local magnetic uniformity criterion
  ! and determines if it is safe to use the drift equations. Estimates new timestep using the Kress [2008]
  ! heuristic or by requiring proper resolution of particle gyromotion. Global parameters drift_epsilon and
  ! orbit_epsilon are used to scale the timesteps to ensure accuracy.
  implicit none

  type(particle) :: myParticle
  real(dp), intent(in) :: tnext
  real(dp), parameter :: kappa_orbit = 1.0d-2 ! These two values are set so that a particle doesn't just
  real(dp), parameter :: kappa_drift = 0.8d-2 ! oscillate between drift and orbit in regions of marginal smoothness
  real(dp), parameter :: dt_base = 1.0d-3
  real(dp) :: b0, wc, kappa, t, igam, t_int, rho
  real(dp), dimension(3) :: bhat, bvec, rvec, rhat
  real(dp), dimension(3) :: fE, fG, fC, fM
  real(dp), dimension(3,3) :: Bgrad
  real(dp), dimension(3) :: gradb, gradbhat, evec, curlb, curlbhat
  real(dp), dimension(:), allocatable :: y
  real(dp), dimension(:), allocatable :: yp
  integer :: i, iflag, nstep, nstep_max

  ! Initialize storage for the integrators and communicators if using RKSuite
  if(myParticle%drift) then
    allocate(y(4),yp(4))
    y = myParticle%y(1:4)
    if(istep==2) call setup(myParticle%comm,myParticle%t,y,tnext,tol,(/(thresh,i=1,4)/),method='H',task='R',message=.FALSE.)
  else
    allocate(y(6),yp(6))
    y = myParticle%y(1:6)
    if(istep==2) call setup(myParticle%comm,myParticle%t,y,tnext,tol,(/(thresh,i=1,6)/),method='M',task='R',message=.FALSE.)
  endif

  nstep = 0
  nstep_max = floor(abs(tnext-myParticle%t)/dt_base)

  do nstep=1,nstep_max

    t = myParticle%t

    ! Check if the timestep is going to take us beyond requested point. If so, set limit of integration appropriately.
    if(itrace>0) then
      t_int = min(tnext,myParticle%t+myParticle%dt)
    else
      t_int = max(tnext,myParticle%t+myParticle%dt)
    endif

    ! Determine the set of equations to solve based on the current drift state of the particle
    if(myParticle%drift) then
      ! Solve guiding center drift equations

      select case(istep)
        case(1)
          ! Use fixed timestep RK4 integrator
          call rk4(drift_derivs,y,t,t_int-t)
          myParticle%y(1:4) = y

        case(2)
          ! Use adaptive RKSuite integrator
          call range_integrate(myParticle%comm,drift_derivs,t_int,t,y,yp,flag=iflag)
          myParticle%t = t
          myParticle%y(1:4) = y
          if(iflag==5) then
            myParticle%integrate = .FALSE.
            write(*,*) "Drift integrator encountered a problem"
          endif

        case default
          call assert(.FALSE.,"stepper_push","option "//int2str(istep)//" not currently supported.")
      end select

    else ! Particle is in full orbit

      select case(istep)
        case(1)
          ! Use fixed timestep RK4 integrator
          call rk4(orbit_derivs,y,t,t_int-t)
          myParticle%y(1:6) = y
        case(2)
          ! Use adaptive RKSuite integrator
          call range_integrate(myParticle%comm,orbit_derivs,t_int,t,y,yp,flag=iflag)
          myParticle%t = t
          myParticle%y(1:6) = y
          if(iflag==5) then
            myParticle%integrate = .FALSE.
            write(*,*) "Orbit integrator encountered a problem"
          endif
        case default
          call assert(.FALSE.,"stepper_push","option "//int2str(istep)//" not currently supported.")
      end select

    endif

    if(.not. myParticle%integrate) exit

    ! Check validity of drift representation
    call get_fields(myParticle,bvec,gradb=Bgrad)

    ! Evaluate the adiabatic drift criterion
    b0 = norm(bvec)
    wc = abs(myParticle%p(1)*b0)
    bhat = bvec/b0
    gradb = matmul(Bgrad,bhat)

    ! Calculate the gyroradius
    if(myParticle%drift) then
      rho = sqrt(2*b0*myParticle%v(2))/wc
    else
      rho = norm(cross_product(myParticle%v,bhat))/wc
    endif

    select case(iswitch)
      case(:-1) ! Guiding center drifts only (GC)
          kappa = 0.9d0*kappa_drift
      case(1:) ! Full particle orbits only (FO)
          kappa = 1.1d0*kappa_drift
      case(0 ) ! Dynamically switch btwn GC and FO
          kappa = rho*norm(gradb)/b0
    end select

    if(kappa > kappa_orbit) then
      ! Use full particle orbits

      if(myParticle%drift) then
        ! Switch to full orbit representation if currently in guiding center
        call drift_to_full(myParticle)

        deallocate(y,yp)
        allocate(y(6),yp(6))

        y = myParticle%y(1:6)

        if(istep==2) then
          ! Modify RKSuite communicator settings for orbit equations
          call collect_garbage(myParticle%comm)
          call setup(myParticle%comm,t,y,tnext,tol,(/(thresh,i=1,6)/),method='M',task='R',message=.FALSE.)
        endif

      endif
    else
      ! Use guiding center equations

      if(.not. myParticle%drift .and. kappa < kappa_drift) then
        ! Switch to guiding center if currently in full particle representation

        call full_to_drift(myParticle)

        deallocate(y,yp)
        allocate(y(4),yp(4))

        y = myParticle%y(1:4)

        if(istep==2) then
          ! Modify RKSuite communicator settings for drift equations
          call collect_garbage(myParticle%comm)
          call setup(myParticle%comm,t,y,tnext,tol,(/(thresh,i=1,4)/),method='H',task='R',message=.FALSE.)
        endif

      endif
    endif

    ! Determine new timestep
    call get_fields(myParticle,bvec,evec=evec,gradb=Bgrad)

    if(myParticle%drift) then
      ! Use Kress-Hudson heuristic
      curlb(1) = Bgrad(2,3)-Bgrad(3,2)
      curlb(2) = Bgrad(3,1)-Bgrad(1,3)
      curlb(3) = Bgrad(1,2)-Bgrad(2,1)
      curlbhat = (curlb-cross_product(gradb,bhat))/b0
      gradbhat = cross_product(curlbhat,bhat)
      igam = 1.d0/sqrt(1.d0+(myParticle%v(1)**2+2.0d0*myParticle%v(2)*b0)/csq)
      fE = -myParticle%p(1)*evec                                              ! Electric force
      fG =  myParticle%v(2)*gradb*igam                                        ! Gradient force
      fC =  myParticle%v(1)**2*gradbhat*igam                                  ! Curvature force
      fM =  myParticle%p(1)*evec-myParticle%v(2)*dot_product(bhat,gradb)*igam ! Mirror force
      myParticle%dt = sign(epsilon_drift*sqrt(myParticle%v(1)**2+2.0d0*myParticle%v(2)*b0)/norm(fE+fG+fC+fM),real(itrace,dp))
    else
      ! Resolve gyro-orbits
      myParticle%dt = sign(epsilon_orbit/wc,real(itrace,dp))
    endif

    ! Truncate timestep if necessary in order to match desired output times
    if(myParticle%t /= tnext) then
      if(itrace<0) then
        if(myParticle%t+myParticle%dt < tnext) myParticle%dt = tnext-myParticle%t
      else
        if(myParticle%t+myParticle%dt > tnext) myParticle%dt = tnext-myParticle%t
      endif
    endif

    ! Check for particle precipitation to ionosphere
    if(norm(myParticle%x) <= 1.0d0) then
      write(*,*) "Particle has precipitated", myParticle%t
      myParticle%integrate=.FALSE.
      exit
    endif

    ! Check if we're prepared to output more data
    if(itrace>0) then
      if(myParticle%t+myParticle%dt > tnext) exit
    else
      if(myParticle%t+myParticle%dt < tnext) exit
    endif

    ! If tracing particle fluxes, terminate integration when it reaches specified boundary
    if(fluxMap .and. (myParticle%x(1) <= xsource)) then
      write(*,*) "Particle has reached source region boundary"
      myParticle%integrate = .FALSE.
      exit
    endif

  enddo

  if(nstep >= nstep_max) then
    write(*,*) "Particle required too many steps for integration"
    myParticle%integrate = .FALSE.
  endif

  ! Clean up integration
  if(istep==2) call collect_garbage(myParticle%comm)
  deallocate(y,yp)

  contains
    ! Making these functions local to the procedure allows us to pass them to the integrator
    ! while maintaining access to the particle information which is necessary to use the
    ! external ODE library RKSuite90

    function orbit_derivs(t,y) result(ydot)
    ! Calculates the derivatives for advancing particle trajectory
    implicit none

    real(dp), dimension(:), intent(in) :: y
    real(dp), dimension(size(y)) :: ydot
    real(dp), intent(in) :: t
    real(dp), dimension(3) :: evec, bvec
    real(dp) :: gami

    if(.not. myParticle%integrate) then
      ! This particle is being taken offline, don't worry about calculating anything.
      ydot = 0.d0
    else

      myParticle%t = t
      myParticle%y(1:6) = y

      call get_fields(myParticle,bvec,evec=evec)
      gami = 1.d0/sqrt(1.d0+dot_product(myParticle%v,myParticle%v)/csq)

      ydot(1:3) =  gami*myParticle%v/re
      ydot(4:6) =  myParticle%p(1)*(evec+gami*cross_product(myParticle%v,bvec))

      if(ndim==2) then
        ! Set z-derivatives to zero
        ydot(3) = 0.d0
        ydot(6) = 0.d0
      endif

    endif

    return

    end function orbit_derivs

    !

    function drift_derivs(t,y) result(ydot)
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

    if(.not. myParticle%integrate) then
      ! This particle has finished integration, take it offline
      ydot = 0.d0
    else
      myParticle%t = t
      myParticle%y(1:4) = y
      call get_fields(myParticle,bvec,evec=evec,gradb=Bgrad,grade=Egrad)
      b0 = norm(bvec)
      bhat = bvec/b0
      gami = 1.d0/sqrt(1.0d0+(2.0d0*b0*myParticle%v(2)+myParticle%v(1)**2)/csq)

      ! Calculate gradient and curl of magnetic field
      gradb = matmul(Bgrad,bhat)

      curlb(1) = Bgrad(2,3)-Bgrad(3,2)
      curlb(2) = Bgrad(3,1)-Bgrad(1,3)
      curlb(3) = Bgrad(1,2)-Bgrad(2,1)

      curle(1) = Egrad(2,3)-Egrad(3,2)
      curle(2) = Egrad(3,1)-Egrad(1,3)
      curle(3) = Egrad(1,2)-Egrad(2,1)

      curlbhat = (curlb-cross_product(gradb,bhat))/b0
      dbhatdt = (bhat*dot_product(bhat,curle)-curle)/b0

      ! Calculate effective fields
      bstar = bvec+myParticle%v(1)*curlbhat/myParticle%p(1)
      bstarpara = dot_product(bhat,bstar)
      estar = evec-(myParticle%v(2)*gradb*gami-myParticle%v(1)*dbhatdt)/myParticle%p(1)

      vdrift = (gami*myParticle%v(1)*bstar+cross_product(estar,bhat))/bstarpara
      dupara = myParticle%p(1)*dot_product(estar,bstar)/bstarpara

      ydot(1:3) = vdrift/re ! Scale drift velocity to re/s
      ydot(4)   = dupara

      if(ndim==2) then ! Set z-derivatives to zero
        ydot(3) = 0.d0
        ydot(4) = 0.d0
      endif

    endif

    return

    end function drift_derivs

  end subroutine stepper_push

!

  subroutine full_to_drift(myParticle)
  ! Switch from full orbit (vx, vy, vz) to guiding center (mu, vpara)

  implicit none

  type(particle) :: myParticle
  real(dp) :: b0, vpara, gam
  real(dp), dimension(3) :: bvec, bhat, vperp

  call get_fields(myParticle,bvec)

  ! Determine local magnetic field and it direction
  b0 = norm(bvec)
  bhat = bvec/b0

  vpara= dot_product(bhat,myParticle%v)
  vperp = myParticle%v-vpara*bhat
  gam = sqrt(1.0+(vpara*vpara+dot_product(vperp,vperp))/csq)

  ! Move the particle to the guiding center
  myParticle%x = myParticle%x + cross_product(myParticle%v,bhat)/(re*b0*myParticle%p(1))

  myParticle%v(1) = vpara
  myParticle%v(2) = dot_product(vperp,vperp)/(2.0d0*b0)

  if(ndim==2) then
    ! Prevent "accidental" introduction of paralell motion or displacements
    myParticle%x(3) = 0.5d0
    myParticle%v(1) = 0.0d0
  endif

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
  real(dp) :: vperp, vpara, rho, wc, b0, phi, dbmin, dbnow
  integer :: i,imin,imax

  call get_fields(myParticle,bvec)
  b0 = norm(bvec)
  bhat = bvec/b0
  wc = abs(b0*myParticle%p(1))
  vpara = myParticle%v(1)
  vperp = sqrt(2.0d0*b0*myParticle%v(2))
  rho = vperp/(re*wc)

  select case(iphase)

    case(1)
      ! Random gyrophase
      R = fac_rotation_matrix(myParticle)
      phi = twopi*random_uniform()

      ! Define unit gyrovector in field-aligned coordinates
      rhat = [cos(phi),sin(phi),0.0d0]

      ! Rotate gyrovector into Cartesian coordinates
      rhat = matmul(R,rhat)
      vhat = sign(1.0d0,myParticle%p(1))*cross_product(bhat,rhat)
      vhat = vhat/norm(vhat)

    case(2)
      ! Brute-force minimization of delta B

      R = fac_rotation_matrix(myParticle,toCartesian=.TRUE.)
      xgc = myParticle%x

      do i=1,nphase
        phi = twopi*real(i,dp)
        rhat = [cos(phi),sin(phi),0.0d0]
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
      vhat = sign(1.0d0,myParticle%p(1))*cross_product(rhat,bhat)
      vhat = vhat/norm(vhat)
      myParticle%x = xgc

  case(3)
    ! Apply Pfefferle's grad B method to minimize delta B

    ! Calculate the local magnetic field and the direction of its gradient
    call get_fields(myParticle,bvec,gradb=Bgrad)

    gradb = matmul(Bgrad,bhat)

    ! Put particle along direction orthogonal to both gradb and b
    rhat = cross_product(bhat,gradb)
    rhat = rhat/norm(rhat)

    vhat = (gradb-bhat*dot_product(bhat,gradb))/norm(cross_product(bvec,gradb))

!    vhat = sign(1.0d0,myParticle%p(1))*cross_product(rhat,bhat)
!    vhat = vhat/norm(vhat)

  end select

  myParticle%x = myParticle%x+rho*rhat
  myParticle%v = vpara*bhat+vperp*vhat

  if(ndim==2) then
    ! Prevent accidental introduction of paralell motion or displacements
    myParticle%x(3) = 0.5d0
    myParticle%v(3) = 0.0d0
  endif

  myParticle%drift = .FALSE.

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

  y = y+(k1+2*k2+2*k3+k4)/6

  return

  end subroutine rk4

end module stepper
