!+MODULE NAME
!   fields
!
!+MODULE DESCRIPTION
!   This module contains the electric and magnetic field arrays and routines that allow us
!   to determine the electric and magnetic fields at arbitrary points in space via interpolation:
!   linear in time and trilinear/tricubic in space
!
!+AUTHOR
!   Jesse Woodroffe
!   jwoodroffe@lanl.gov

module fields

use global
use finite_differences
use interpolation

contains

  subroutine fields_initialize()
  ! Read in and initialize the global E and B profiles
  implicit none

  integer :: i,j
  character(len=22) :: dfile

  ! Allocate storage for global arrays
  allocate(xgrid(nx),ygrid(ny))

  if(ndim==2) then
    allocate(zgrid(2))
    allocate(BX3D(nt,nx,ny,2),BY3D(nt,nx,ny,2),BZ3D(nt,nx,ny,2))
    allocate(EX3D(nt,nx,ny,2),EY3D(nt,nx,ny,2),EZ3D(nt,nx,ny,2))
  else
    allocate(zgrid(nz))
    allocate(BX3D(nt,nx,ny,nz),BY3D(nt,nx,ny,nz),BZ3D(nt,nx,ny,nz))
    allocate(EX3D(nt,nx,ny,nz),EY3D(nt,nx,ny,nz),EZ3D(nt,nx,ny,nz))
  endif

  ! Read data files files
  call read_array('ptm_data/xgrid.bin',xgrid)
  call read_array('ptm_data/ygrid.bin',ygrid)

  if(ndim==3) then
    call read_array('ptm_data/zgrid.bin',zgrid)
  else
    zgrid = (/0.d0,1.d0/)
  endif

  ! Spatial and temporal ranges
  if(dtIn > 0.0d0) then
    ! We are specifying what files to use, determine the times corresponding to those files

    allocate(tgrid(nt))
    TMin = 0.0d0
    TMax = dtIn*real(nt-1,dp)
    tgrid = linspace(TMin,Tmax,nt)

  else

    ! We are reading in the grid directly, figure out what files are needed

    allocate(tgrid(ntot))

    call read_array('ptm_data/tgrid.bin',tgrid)

    ifirst = maxval(maxloc(tgrid,tgrid<=tlo))
    ilast = maxval(minloc(tgrid,tgrid>=thi))

    nt = ilast-ifirst+1

    TMin = tgrid(ifirst)
    TMax = tgrid(ilast)

  endif

  XMin = minval(xgrid)
  XMax = maxval(xgrid)
  YMin = minval(ygrid)
  YMax = maxval(ygrid)
  ZMin = minval(zgrid)
  ZMax = maxval(zgrid)

  write(*,*) "SIMULATION DOMAIN"
  write(*,*) xmin, " <= X <= ", xmax
  write(*,*) ymin, " <= Y <= ", ymax
  write(*,*) zmin, " <= Z <= ", zmax

  ! Get electromagnetic fields

  if(ndim==2) then

    write(*,*) "Reading 2D electromagnetic field data"

    do j=ifirst,ilast
      i=j-ifirst+1

      write(dfile,'("ptm_data/bx2d_",i4.4,".bin")') j
      call read_array(dfile,BX3D(i,:,:,1))
      BX3D(i,:,:,2) = BX3D(i,:,:,1)

      write(dfile,'("ptm_data/by2d_",i4.4,".bin")') j
      call read_array(dfile,BY3D(i,:,:,1))
      BY3D(i,:,:,2) = BY3D(i,:,:,1)

      write(dfile,'("ptm_data/bz2d_",i4.4,".bin")') j
      call read_array(dfile,BZ3D(i,:,:,1))
      BZ3D(i,:,:,2) = BZ3D(i,:,:,1)

      write(dfile,'("ptm_data/ex2d_",i4.4,".bin")') j
      call read_array(dfile,EX3D(i,:,:,1))
      EX3D(i,:,:,2) = EX3D(i,:,:,1)

      write(dfile,'("ptm_data/ey2d_",i4.4,".bin")') j
      call read_array(dfile,EY3D(i,:,:,1))
      EY3D(i,:,:,2) = EY3D(i,:,:,1)

      write(dfile,'("ptm_data/ez2d_",i4.4,".bin")') j
      call read_array(dfile,EZ3D(i,:,:,1))
      EZ3D(i,:,:,2) = EZ3D(i,:,:,1)

    enddo

  else

    write(*,*) "Reading 3D electromagnetic field data"

    do j=ifirst,ilast

      i = j-ifirst+1

      write(dfile,'("ptm_data/bx3d_",i4.4,".bin")') j

      call read_array(dfile,BX3D(i,:,:,:))

      write(dfile,'("ptm_data/by3d_",i4.4,".bin")') j
      call read_array(dfile,BY3D(i,:,:,:))

      write(dfile,'("ptm_data/bz3d_",i4.4,".bin")') j
      call read_array(dfile,BZ3D(i,:,:,:))

      write(dfile,'("ptm_data/ex3d_",i4.4,".bin")') j
      call read_array(dfile,EX3D(i,:,:,:))

      write(dfile,'("ptm_data/ey3d_",i4.4,".bin")') j
      call read_array(dfile,EY3D(i,:,:,:))

      write(dfile,'("ptm_data/ez3d_",i4.4,".bin")') j
      call read_array(dfile,EZ3D(i,:,:,:))

    enddo

  endif

  return

  end subroutine fields_initialize

!

  subroutine get_fields(myParticle,bvec,evec,gradb,grade,doInit)
  ! Calculate the magnetic field vector at a point and, optionally, the
  ! electric field vector and gradients of each magnetic field component
  implicit none
  type(particle) :: myParticle
  real(dp), dimension(3), intent(out) :: bvec
  real(dp), dimension(3), intent(out), optional :: evec
  real(dp), dimension(3,3), intent(out), optional :: gradb, grade
  real(dp), dimension(3) ::  grad1, grad2, xseg, yseg, zseg, bhat
  real(dp) :: dx, dy, dz, dt
  real(dp), dimension(2,2,2) :: bx_x, bx_y, bx_z, bx_xy, bx_xz, bx_yz, bx_xyz
  real(dp), dimension(2,2,2) :: by_x, by_y, by_z, by_xy, by_xz, by_yz, by_xyz
  real(dp), dimension(2,2,2) :: bz_x, bz_y, bz_z, bz_xy, bz_xz, bz_yz, bz_xyz
  real(dp), dimension(2,2,2) :: ex_x, ex_y, ex_z, ex_xy, ex_xz, ex_yz, ex_xyz
  real(dp), dimension(2,2,2) :: ey_x, ey_y, ey_z, ey_xy, ey_xz, ey_yz, ey_xyz
  real(dp), dimension(2,2,2) :: ez_x, ez_y, ez_z, ez_xy, ez_xz, ez_yz, ez_xyz
  real(dp) :: b1, b2, e1, e2
  real(dp) :: xi, xim
  integer :: im, jm, km, lm, l, i, j, k, ii, jj, kk, ll
  integer :: imin, imax, jmin, jmax, kmin, kmax
  logical, intent(in), optional :: doInit
  integer :: idx, idy, idz

  ! Check if particle is still in the 4D simulation domain
  if(any((/myParticle%t    <  tmin,myParticle%t     > tmax, &
           myParticle%x(1) <= xmin,myParticle%x(1) >= xmax, &
           myParticle%x(2) <= ymin,myParticle%x(2) >= ymax, &
           myParticle%x(3) <= zmin,myParticle%x(3) >= zmax/))) then

    if(myParticle%x(1) < xmin .or. myParticle%x(1) > xmax) write(*,'(a20,3es16.5)') "X OUT OF BOUNDS", xmin, myParticle%x(1), xmax
    if(myParticle%x(2) < ymin .or. myParticle%x(2) > ymax) write(*,'(a20,3es16.5)') "Y OUT OF BOUNDS", ymin, myParticle%x(2), ymax
    if(myParticle%x(3) < zmin .or. myParticle%x(3) > zmax) write(*,'(a20,3es16.5)') "Z OUT OF BOUNDS", zmin, myParticle%x(3), zmax

    write(*,*) "PARTICLE POSITION WHEN LEAVING DOMAIN"
    write(*,*) tmin, myparticle%t,    tmax
    write(*,*) xmin, myparticle%x(1), xmax
    write(*,*) ymin, myparticle%x(2), ymax
    write(*,*) zmin, myparticle%x(3), zmax

    ! Mark particle as finished
    myParticle%integrate=.FALSE.

    ! Put the particle on the nearest edge of the 4D domain for flux mapping
    ! (These values won't be used for further trajectory integration)

    ! ::NOTE THE COMMENTED SECTIONS ON THE FOLLOWING LINES, NEED TO CONFIRM THIS DOESN'T BREAK ANYTHING::

    if(myParticle%x(1) <= xmin) myParticle%x(1) = xmin!+0.5*myParticle%grid%dx
    if(myParticle%x(1) >= xmax) myParticle%x(1) = xmax!-0.5*myParticle%grid%dx
    if(myParticle%x(2) <= ymin) myParticle%x(2) = ymin!+0.5*myParticle%grid%dy
    if(myParticle%x(2) >= ymax) myParticle%x(2) = ymax!-0.5*myParticle%grid%dy
    if(ndim==2) then
      myParticle%x(3) = 0.5d0
    else
      if(myParticle%x(3) <= zmin) myParticle%x(3) = zmin+0.5*myParticle%grid%dz
      if(myParticle%x(3) >= zmax) myParticle%x(3) = zmax-0.5*myParticle%grid%dz
    endif
    if(myParticle%t < tmin) myParticle%t = tmin
    if(myParticle%t > tmax) myParticle%t = tmax

  endif

  im = locate(xgrid,myParticle%x(1))
  jm = locate(ygrid,myParticle%x(2))
  km = locate(zgrid,myParticle%x(3))
!  lm = locate(tgrid,myParticle%t)

  if(TMin < minval(tgrid) .and. myParticle%t < minval(tgrid)) then
    ! We assume fields remain constant at earlier times
    lm = 1
  else if (TMax > maxval(tgrid) .and. myParticle%t > maxval(tgrid)) then
    ! We assume fields remain constant at later times
    lm = nt-1
  else
    lm = locate(tgrid,myParticle%t)
  endif

  ! Prevent segfaults when particles are exactly on simulation boundaries
  ! (Mainly an issue with time during backwards tracing, dt calculation causes segfault if lm is not decremented here)
  if(im==nx) im = im-1
  if(jm==ny) jm = jm-1
  if(km==nz) km = km-1
  if(lm==nt) lm = lm-1

  ! Update the grid and initialize the interpolators
  if(any((/im .ne. myParticle%grid%im, &
           jm .ne. myParticle%grid%jm, &
           km .ne. myParticle%grid%km, &
           lm .ne. myParticle%grid%lm, &
           present(doInit)/))) then

    ! Update grid spacings
    dx = xgrid(im+1)-xgrid(im)
    dy = ygrid(jm+1)-ygrid(jm)
    dz = zgrid(km+1)-zgrid(km)
    dt = tgrid(lm+1)-tgrid(lm)

    ! Configure the 4d grid
    call grid_init(myParticle%grid,im,jm,km,lm,dx,dy,dz,dt)

    ! Approximate derivatives at grid corners for tricubic interpolator
    do ll=1,2
      l=lm+ll-1
      do ii=1,2
        i=im+ii-1
        if(i==1) then
          idx = -1
          imin = 1
          imax = 3
        else if(i==nx) then
          idx = 1
          imin = nx-2
          imax = nx
        else
          idx = 0
          imin = i-1
          imax = i+1
        endif
        xseg=xgrid(imin:imax)
        do jj=1,2
          j = jm+jj-1
          if(j==1) then
            idy = -1
            jmin = 1
            jmax = 3
          else if(j==ny) then
            idy = 1
            jmin = ny-2
            jmax = ny
          else
            idy = 0
            jmin = j-1
            jmax = j+1
          endif
          yseg=ygrid(jmin:jmax)
          do kk=1,2
            k=km+kk-1
            if(k==1) then
              idz = -1
              kmin = 1
              kmax = 3
            else if(k==nz) then
              idz = 1
              kmin = nz-2
              kmax = nz
            else
              idz = 0
              kmin = k-1
              kmax = k+1
            endif
            zseg=zgrid(kmin:kmax)

            !**********************************
            ! DERIVATIVES OF THE MAGNETIC FIELD

            ! Derivatives of the X-component of the field
            bx_x(ii,jj,kk)  = deriv1(xseg,bx3d(l,imin:imax,j,k),idx=idx)
            bx_y(ii,jj,kk)  = deriv1(yseg,bx3d(l,i,jmin:jmax,k),idx=idy)
            bx_xy(ii,jj,kk) = deriv2(xseg,yseg,bx3d(l,imin:imax,jmin:jmax,k),idx=idx,idy=idy)

            if(ndim==2) then ! Uniform in the z-direction
              bx_z(ii,jj,kk)   = 0.d0
              bx_xz(ii,jj,kk)  = 0.d0
              bx_yz(ii,jj,kk)  = 0.d0
              bx_xyz(ii,jj,kk) = 0.d0
            else
              bx_z(ii,jj,kk)   = deriv1(zseg,bx3d(l,i,j,kmin:kmax),idx=idz)
              bx_xz(ii,jj,kk)  = deriv2(xseg,zseg,bx3d(l,imin:imax,j,kmin:kmax),idx=idx,idy=idz)
              bx_yz(ii,jj,kk)  = deriv2(yseg,zseg,bx3d(l,i,jmin:jmax,kmin:kmax),idx=idy,idy=idz)
              bx_xyz(ii,jj,kk) = deriv3(xseg,yseg,zseg,bx3d(l,imin:imax,jmin:jmax,kmin:kmax),idx=idx,idy=idy,idz=idz)
            endif

            ! Derivatives of the X-component of the field
            by_x(ii,jj,kk)  = deriv1(xseg,by3d(l,imin:imax,j,k),idx=idx)
            by_y(ii,jj,kk)  = deriv1(yseg,by3d(l,i,jmin:jmax,k),idx=idy)
            by_xy(ii,jj,kk) = deriv2(xseg,yseg,by3d(l,imin:imax,jmin:jmax,k),idx=idx,idy=idy)

            if(ndim==2) then ! Uniform in the z-direction
              by_z(ii,jj,kk)   = 0.d0
              by_xz(ii,jj,kk)  = 0.d0
              by_yz(ii,jj,kk)  = 0.d0
              by_xyz(ii,jj,kk) = 0.d0
            else
              by_z(ii,jj,kk)   = deriv1(zseg,by3d(l,i,j,kmin:kmax),idx=idz)
              by_xz(ii,jj,kk)  = deriv2(xseg,zseg,by3d(l,imin:imax,j,kmin:kmax),idx=idx,idy=idz)
              by_yz(ii,jj,kk)  = deriv2(yseg,zseg,by3d(l,i,jmin:jmax,kmin:kmax),idx=idy,idy=idz)
              by_xyz(ii,jj,kk) = deriv3(xseg,yseg,zseg,by3d(l,imin:imax,jmin:jmax,kmin:kmax),idx=idx,idy=idy,idz=idz)
            endif

            ! Derivatives of the X-component of the field
            bz_x(ii,jj,kk)  = deriv1(xseg,bz3d(l,imin:imax,j,k),idx=idx)
            bz_y(ii,jj,kk)  = deriv1(yseg,bz3d(l,i,jmin:jmax,k),idx=idy)
            bz_xy(ii,jj,kk) = deriv2(xseg,yseg,bz3d(l,imin:imax,jmin:jmax,k),idx=idx,idy=idy)

            if(ndim==2) then ! Uniform in the z-direction
              bz_z(ii,jj,kk)   = 0.d0
              bz_xz(ii,jj,kk)  = 0.d0
              bz_yz(ii,jj,kk)  = 0.d0
              bz_xyz(ii,jj,kk) = 0.d0
            else
              bz_z(ii,jj,kk)   = deriv1(zseg,bz3d(l,i,j,kmin:kmax),idx=idz)
              bz_xz(ii,jj,kk)  = deriv2(xseg,zseg,bz3d(l,imin:imax,j,kmin:kmax),idx=idx,idy=idz)
              bz_yz(ii,jj,kk)  = deriv2(yseg,zseg,bz3d(l,i,jmin:jmax,kmin:kmax),idx=idy,idy=idz)
              bz_xyz(ii,jj,kk) = deriv3(xseg,yseg,zseg,bz3d(l,imin:imax,jmin:jmax,kmin:kmax),idx=idx,idy=idy,idz=idz)
            endif

            !**********************************
            ! DERIVATIVES OF THE ELECTRIC FIELD

            ! Derivatives of the X-component
            ex_x(ii,jj,kk)  = deriv1(xseg,ex3d(l,imin:imax,j,k),idx=idx)
            ex_y(ii,jj,kk)  = deriv1(yseg,ex3d(l,i,jmin:jmax,k),idx=idy)
            ex_xy(ii,jj,kk) = deriv2(xseg,yseg,ex3d(l,imin:imax,jmin:jmax,k),idx=idx,idy=idy)

            if(ndim==2) then ! Uniform in the z-direction
              ex_z(ii,jj,kk)   = 0.d0
              ex_xz(ii,jj,kk)  = 0.d0
              ex_yz(ii,jj,kk)  = 0.d0
              ex_xyz(ii,jj,kk) = 0.d0
            else
              ex_z(ii,jj,kk)   = deriv1(zseg,ex3d(l,i,j,kmin:kmax),idx=idz)
              ex_xz(ii,jj,kk)  = deriv2(xseg,zseg,ex3d(l,imin:imax,j,kmin:kmax),idx=idx,idy=idz)
              ex_yz(ii,jj,kk)  = deriv2(yseg,zseg,ex3d(l,i,jmin:jmax,kmin:kmax),idx=idy,idy=idz)
              ex_xyz(ii,jj,kk) = deriv3(xseg,yseg,zseg,ex3d(l,imin:imax,jmin:jmax,kmin:kmax),idx=idx,idy=idy,idz=idz)
            endif

            ! Derivatives of the Y-component
            ey_x(ii,jj,kk)  = deriv1(xseg,ey3d(l,imin:imax,j,k),idx=idx)
            ey_y(ii,jj,kk)  = deriv1(yseg,ey3d(l,i,jmin:jmax,k),idx=idy)
            ey_xy(ii,jj,kk) = deriv2(xseg,yseg,ey3d(l,imin:imax,jmin:jmax,k),idx=idx,idy=idy)

            if(ndim==2) then ! Uniform in the z-direction
              ey_z(ii,jj,kk)   = 0.d0
              ey_xz(ii,jj,kk)  = 0.d0
              ey_yz(ii,jj,kk)  = 0.d0
              ey_xyz(ii,jj,kk) = 0.d0
            else
              ey_z(ii,jj,kk)   = deriv1(zseg,ey3d(l,i,j,kmin:kmax),idx=idz)
              ey_xz(ii,jj,kk)  = deriv2(xseg,zseg,ey3d(l,imin:imax,j,kmin:kmax),idx=idx,idy=idz)
              ey_yz(ii,jj,kk)  = deriv2(yseg,zseg,ey3d(l,i,jmin:jmax,kmin:kmax),idx=idy,idy=idz)
              ey_xyz(ii,jj,kk) = deriv3(xseg,yseg,zseg,ey3d(l,imin:imax,jmin:jmax,kmin:kmax),idx=idx,idy=idy,idz=idz)
            endif

            ! Derivatives of the Z-component
            ez_x(ii,jj,kk)  = deriv1(xseg,ez3d(l,imin:imax,j,k),idx=idx)
            ez_y(ii,jj,kk)  = deriv1(yseg,ez3d(l,i,jmin:jmax,k),idx=idy)
            ez_xy(ii,jj,kk) = deriv2(xseg,yseg,ez3d(l,imin:imax,jmin:jmax,k),idx=idx,idy=idy)

            if(ndim==2) then ! Uniform in the z-direction
              ez_z(ii,jj,kk)   = 0.d0
              ez_xz(ii,jj,kk)  = 0.d0
              ez_yz(ii,jj,kk)  = 0.d0
              ez_xyz(ii,jj,kk) = 0.d0
            else
              ez_z(ii,jj,kk)   = deriv1(zseg,ez3d(l,i,j,kmin:kmax),idx=idz)
              ez_xz(ii,jj,kk)  = deriv2(xseg,zseg,ez3d(l,imin:imax,j,kmin:kmax),idx=idx,idy=idz)
              ez_yz(ii,jj,kk)  = deriv2(yseg,zseg,ez3d(l,i,jmin:jmax,kmin:kmax),idx=idy,idy=idz)
              ez_xyz(ii,jj,kk) = deriv3(xseg,yseg,zseg,ez3d(l,imin:imax,jmin:jmax,kmin:kmax),idx=idx,idy=idy,idz=idz)
            endif

          enddo
        enddo
      enddo

      ! Get the tricubic interpolation coefficients for all requested fields
      call tricubic_init(myParticle%bxinterp(ll,:), myParticle%grid, BX3D(l,im:im+1,jm:jm+1,km:km+1),BX_X,BX_Y,BX_Z,BX_XY,BX_XZ,BX_YZ,BX_XYZ)
      call tricubic_init(myParticle%byinterp(ll,:), myParticle%grid, BY3D(l,im:im+1,jm:jm+1,km:km+1),BY_X,BY_Y,BY_Z,BY_XY,BY_XZ,BY_YZ,BY_XYZ)
      call tricubic_init(myParticle%bzinterp(ll,:), myParticle%grid, BZ3D(l,im:im+1,jm:jm+1,km:km+1),BZ_X,BZ_Y,BZ_Z,BZ_XY,BZ_XZ,BZ_YZ,BZ_XYZ)

      call tricubic_init(myParticle%exinterp(ll,:), myParticle%grid, EX3D(l,im:im+1,jm:jm+1,km:km+1),EX_X,EX_Y,EX_Z,EX_XY,EX_XZ,EX_YZ,EX_XYZ)
      call tricubic_init(myParticle%eyinterp(ll,:), myParticle%grid, EY3D(l,im:im+1,jm:jm+1,km:km+1),EY_X,EY_Y,EY_Z,EY_XY,EY_XZ,EY_YZ,EY_XYZ)
      call tricubic_init(myParticle%ezinterp(ll,:), myParticle%grid, EZ3D(l,im:im+1,jm:jm+1,km:km+1),EZ_X,EZ_Y,EZ_Z,EZ_XY,EZ_XZ,EZ_YZ,EZ_XYZ)

    enddo
  endif

  !**********************************************
  ! Coefficients for linear interpolation in time

  xi = (myParticle%t-tgrid(myParticle%grid%lm))/myParticle%grid%dt
  xim = 1.0-xi

  ! Interpolation and gradient calculation for the electric field
  if(present(evec)) then
    if(present(grade)) then
      call tricubic_interpolate(myParticle%exinterp(1,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),e1,gradf=grad1)
      call tricubic_interpolate(myParticle%exinterp(2,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),e2,gradf=grad2)
      grade(:,1) = (xim*grad1+xi*grad2)/re
      evec(1) = xim*e1+xi*e2
      call tricubic_interpolate(myParticle%eyinterp(1,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),e1,gradf=grad1)
      call tricubic_interpolate(myParticle%eyinterp(2,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),e2,gradf=grad2)
      grade(:,2) = (xim*grad1+xi*grad2)/re
      evec(2) = xim*e1+xi*e2
      call tricubic_interpolate(myParticle%ezinterp(1,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),e1,gradf=grad1)
      call tricubic_interpolate(myParticle%ezinterp(2,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),e2,gradf=grad2)
      grade(:,3) = (xim*grad1+xi*grad2)/re
      evec(3) = xim*e1+xi*e2
    else
      call tricubic_interpolate(myParticle%exinterp(1,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),e1)
      call tricubic_interpolate(myParticle%exinterp(2,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),e2)
      evec(1) = xim*e1+xi*e2
      call tricubic_interpolate(myParticle%eyinterp(1,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),e1)
      call tricubic_interpolate(myParticle%eyinterp(2,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),e2)
      evec(2) = xim*e1+xi*e2
      call tricubic_interpolate(myParticle%ezinterp(1,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),e1)
      call tricubic_interpolate(myParticle%ezinterp(2,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),e2)
      evec(3) = xim*e1+xi*e2
    endif
  endif

  ! Interpolation and gradient calculation for the magnetic field
  if(present(gradb)) then
    call tricubic_interpolate(myParticle%bxinterp(1,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),b1,gradf=grad1)
    call tricubic_interpolate(myParticle%bxinterp(2,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),b2,gradf=grad2)
    gradb(:,1) = (xim*grad1+xi*grad2)/re
    bvec(1) = xim*b1+xi*b2
    call tricubic_interpolate(myParticle%byinterp(1,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),b1,gradf=grad1)
    call tricubic_interpolate(myParticle%byinterp(2,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),b2,gradf=grad2)
    gradb(:,2) = (xim*grad1+xi*grad2)/re
    bvec(2) = xim*b1+xi*b2
    call tricubic_interpolate(myParticle%bzinterp(1,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),b1,gradf=grad1)
    call tricubic_interpolate(myParticle%bzinterp(2,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),b2,gradf=grad2)
    gradb(:,3) = (xim*grad1+xi*grad2)/re
    bvec(3) = xim*b1+xi*b2
  else
    call tricubic_interpolate(myParticle%bxinterp(1,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),b1)
    call tricubic_interpolate(myParticle%bxinterp(2,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),b2)
    bvec(1) = xim*b1+xi*b2
    call tricubic_interpolate(myParticle%byinterp(1,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),b1)
    call tricubic_interpolate(myParticle%byinterp(2,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),b2)
    bvec(2) = xim*b1+xi*b2
    call tricubic_interpolate(myParticle%bzinterp(1,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),b1)
    call tricubic_interpolate(myParticle%bzinterp(2,:),myParticle%grid,myParticle%x(1),myParticle%x(2),myParticle%x(3),b2)
    bvec(3) = xim*b1+xi*b2
  endif

  ! Remove parallel electric fields that were introduced by interpolation [J. Birn]
  if(present(evec) .and. norm(bvec) > thresh) then
    bhat = bvec/norm(bvec)
    evec = evec-dot_product(evec,bhat)*bhat
  endif

  return

  end subroutine get_fields

!

  subroutine fields_cleanup()
  ! Free up allocated storage related to this module
  implicit none

  deallocate(xgrid,ygrid,zgrid,tgrid)
  deallocate(BX3D,BY3D,BZ3D)
  deallocate(EX3D,EY3D,EZ3D)

  return

  end subroutine fields_cleanup

end module fields
