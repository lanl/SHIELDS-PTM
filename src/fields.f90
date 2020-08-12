!+MODULE NAME
!   fields
!
!+MODULE DESCRIPTION
!   This module contains the electric and magnetic field arrays and routines that allow us
!   to determine the electric and magnetic fields at arbitrary points in space via interpolation:
!   linear in time and trilinear/tricubic in space
!
!+CONTAINS subroutines
!   FIELD_INITIALIZE
!   GET_FIELDS
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
  character(len=28) :: dfile

  ! Read temporal and spatial grids

  call read_grid('ptm_data/tgrid.dat',tgrid)
  ifirst = maxval(maxloc(tgrid,tgrid<=tlo))
  ilast = maxval(maxloc(tgrid,tgrid>=thi))
  nt = ilast-ifirst+1
  TMin = tgrid(ifirst)
  TMax = tgrid(ilast)

  write(dfile,'("ptm_data/ptm_fields_",i4.4,".dat")') ifirst
  call read_fields(dfile,xv=xgrid,yv=ygrid,zv=zgrid)

  nx = size(xgrid)
  ny = size(ygrid)
  nz = size(zgrid)

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
  write(*,*) tmin, " <= T <= ", tmax

  ! Allocate storage for electromagnetic fields
  allocate(BX3D(nt,nx,ny,nz))
  allocate(BY3D,BZ3D,EX3D,EY3D,EZ3D,mold=BX3D)

  ! Get electromagnetic fields
  write(*,*) "Reading 3D electromagnetic field data"

  do j=ifirst,ilast

    i = j-ifirst+1

    write(dfile,'("ptm_data/ptm_fields_",i4.4,".dat")') j
    call read_fields(dfile,bxm=BX3D(i,:,:,:),bym=BY3D(i,:,:,:),bzm=BZ3D(i,:,:,:), &
                           exm=EX3D(i,:,:,:),eym=EY3D(i,:,:,:),ezm=EZ3D(i,:,:,:))

  enddo

  return

  end subroutine fields_initialize

!

  subroutine read_fields(fname,xv,yv,zv,bxm,bym,bzm,exm,eym,ezm)
  ! Read data from a regular cartesian grid in PTM's XYZ format. The general structure of the file is
  !
  ! nx ny nz
  ! x(1) x(2) ... x(nx-1) x(nx)
  ! y(1) y(2) ... y(ny-1) y(ny)
  ! z(1) z(2) ... z(nz-1) z(nz)
  ! i1 j1 k1 bx(1) by(1) bz(1) ex(1) ey(1) ez(1)
  ! ...
  ! ...
  ! iN jN kN bx(N) by(N) bz(N) ex(N) ey(N) ez(N)
  !
  ! Where i, j, k are indices corresponding to a particular element of the x, y, or z grid respectively
  ! (i.e., the logical position) and N = nx*ny*nz.
  !
  ! The purpose of this routine is to make the field-specification data files more readable and flexible.
  ! This routine allows you to pick single components (e.g. ey only) or all of them, but the calling
  ! sequence requires each argument other than fname to be specified as a keyword.

  implicit none

  character(len=*), intent(in) :: fname
  real(dp), dimension(:), allocatable, optional :: xv, yv, zv
  real(dp), dimension(:,:,:), optional :: bxm, bym, bzm, exm, eym, ezm
  real(dp) :: bx, by, bz, ex, ey, ez
  integer :: lun, i, j, k, istat, n, ntot, numx, numy, numz

  open(newunit=lun,file=fname,status='old',iostat=istat,action='read')
  call assert(istat==0,'read_fields','Error opening '//fname)

  read(lun,*) numx, numy, numz

  ntot = numx*numy*numz

  ! Allocate and read the grids if appropriate, otherwise just skip past them
  if(present(xv)) then
      if(.not. allocated(xv)) allocate(xv(numx))
      call assert(size(xv)==numx,'read_fields','xgrid size mismatch')
      read(lun,*) xv
  else
      read(lun,*)
  endif

  if(present(yv)) then
      if(.not. allocated(yv)) allocate(yv(numy))
      call assert(size(yv)==numy,'read_fields','ygrid size mismatch')
      read(lun,*) yv
  else
      read(lun,*)
  endif

  if(present(zv)) then
      if(.not. allocated(zv)) allocate(zv(numz))
      call assert(size(zv)==numz,'read_fields','zgrid size mismatch')
      read(lun,*) zv
  else
      read(lun,*)
  endif

  ! If user has asked for any field components, read the file and grab the appropriate ones.
  if(any([present(bxm),present(bym),present(bzm),present(exm),present(eym),present(ezm)])) then

    ! Check that fields are properly allocated and have the appropriate shape to hold the data
    if(present(bxm)) call assert(all([size(bxm,1)==numx,size(bxm,2)==numy,size(bxm,3)==numz]),'read_fields','bx shape mismatch')
    if(present(bym)) call assert(all([size(bym,1)==numx,size(bym,2)==numy,size(bym,3)==numz]),'read_fields','by shape mismatch')
    if(present(bzm)) call assert(all([size(bzm,1)==numx,size(bzm,2)==numy,size(bzm,3)==numz]),'read_fields','bz shape mismatch')
    if(present(exm)) call assert(all([size(exm,1)==numx,size(exm,2)==numy,size(exm,3)==numz]),'read_fields','ex shape mismatch')
    if(present(eym)) call assert(all([size(eym,1)==numx,size(eym,2)==numy,size(eym,3)==numz]),'read_fields','ey shape mismatch')
    if(present(ezm)) call assert(all([size(ezm,1)==numx,size(ezm,2)==numy,size(ezm,3)==numz]),'read_fields','ez shape mismatch')

    do n=1,ntot
      read(lun,*,iostat=istat) i, j, k, bx, by, bz, ex, ey, ez
      call assert(istat==0,'read_fields','Error reading '//fname)

      if(present(bxm)) bxm(i,j,k)=bx
      if(present(bym)) bym(i,j,k)=by
      if(present(bzm)) bzm(i,j,k)=bz
      if(present(exm)) exm(i,j,k)=ex
      if(present(eym)) eym(i,j,k)=ey
      if(present(ezm)) ezm(i,j,k)=ez

    enddo

  endif

  close(lun)

  return

  end subroutine read_fields

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
    write(*,*) 'T:',tmin, myparticle%t,    tmax
    write(*,*) 'X:',xmin, myparticle%x(1), xmax
    write(*,*) 'Y:',ymin, myparticle%x(2), ymax
    write(*,*) 'Z:',zmin, myparticle%x(3), zmax

    ! Mark particle as finished
    myParticle%integrate=.FALSE.

    ! Put the particle on the nearest edge of the 4D domain for flux mapping
    ! (These values won't be used for further trajectory integration)

    ! ::NOTE THE COMMENTED SECTIONS ON THE FOLLOWING LINES, NEED TO CONFIRM THIS DOESN'T BREAK ANYTHING::

    if(myParticle%x(1) <= xmin) myParticle%x(1) = xmin!+0.5*myParticle%grid%dx
    if(myParticle%x(1) >= xmax) myParticle%x(1) = xmax!-0.5*myParticle%grid%dx
    if(myParticle%x(2) <= ymin) myParticle%x(2) = ymin!+0.5*myParticle%grid%dy
    if(myParticle%x(2) >= ymax) myParticle%x(2) = ymax!-0.5*myParticle%grid%dy
    if(myParticle%x(3) <= zmin) myParticle%x(3) = zmin!+0.5*myParticle%grid%dz
    if(myParticle%x(3) >= zmax) myParticle%x(3) = zmax!-0.5*myParticle%grid%dz
    if(myParticle%t < tmin) myParticle%t = tmin
    if(myParticle%t > tmax) myParticle%t = tmax

  endif

  im = locate(xgrid,myParticle%x(1))
  jm = locate(ygrid,myParticle%x(2))
  km = locate(zgrid,myParticle%x(3))

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

            bx_z(ii,jj,kk)   = deriv1(zseg,bx3d(l,i,j,kmin:kmax),idx=idz)
            bx_xz(ii,jj,kk)  = deriv2(xseg,zseg,bx3d(l,imin:imax,j,kmin:kmax),idx=idx,idy=idz)
            bx_yz(ii,jj,kk)  = deriv2(yseg,zseg,bx3d(l,i,jmin:jmax,kmin:kmax),idx=idy,idy=idz)
            bx_xyz(ii,jj,kk) = deriv3(xseg,yseg,zseg,bx3d(l,imin:imax,jmin:jmax,kmin:kmax),idx=idx,idy=idy,idz=idz)

            ! Derivatives of the X-component of the field
            by_x(ii,jj,kk)  = deriv1(xseg,by3d(l,imin:imax,j,k),idx=idx)
            by_y(ii,jj,kk)  = deriv1(yseg,by3d(l,i,jmin:jmax,k),idx=idy)
            by_xy(ii,jj,kk) = deriv2(xseg,yseg,by3d(l,imin:imax,jmin:jmax,k),idx=idx,idy=idy)

            by_z(ii,jj,kk)   = deriv1(zseg,by3d(l,i,j,kmin:kmax),idx=idz)
            by_xz(ii,jj,kk)  = deriv2(xseg,zseg,by3d(l,imin:imax,j,kmin:kmax),idx=idx,idy=idz)
            by_yz(ii,jj,kk)  = deriv2(yseg,zseg,by3d(l,i,jmin:jmax,kmin:kmax),idx=idy,idy=idz)
            by_xyz(ii,jj,kk) = deriv3(xseg,yseg,zseg,by3d(l,imin:imax,jmin:jmax,kmin:kmax),idx=idx,idy=idy,idz=idz)

            ! Derivatives of the X-component of the field
            bz_x(ii,jj,kk)  = deriv1(xseg,bz3d(l,imin:imax,j,k),idx=idx)
            bz_y(ii,jj,kk)  = deriv1(yseg,bz3d(l,i,jmin:jmax,k),idx=idy)
            bz_xy(ii,jj,kk) = deriv2(xseg,yseg,bz3d(l,imin:imax,jmin:jmax,k),idx=idx,idy=idy)

            bz_z(ii,jj,kk)   = deriv1(zseg,bz3d(l,i,j,kmin:kmax),idx=idz)
            bz_xz(ii,jj,kk)  = deriv2(xseg,zseg,bz3d(l,imin:imax,j,kmin:kmax),idx=idx,idy=idz)
            bz_yz(ii,jj,kk)  = deriv2(yseg,zseg,bz3d(l,i,jmin:jmax,kmin:kmax),idx=idy,idy=idz)
            bz_xyz(ii,jj,kk) = deriv3(xseg,yseg,zseg,bz3d(l,imin:imax,jmin:jmax,kmin:kmax),idx=idx,idy=idy,idz=idz)

            !**********************************
            ! DERIVATIVES OF THE ELECTRIC FIELD

            ! Derivatives of the X-component
            ex_x(ii,jj,kk)  = deriv1(xseg,ex3d(l,imin:imax,j,k),idx=idx)
            ex_y(ii,jj,kk)  = deriv1(yseg,ex3d(l,i,jmin:jmax,k),idx=idy)
            ex_xy(ii,jj,kk) = deriv2(xseg,yseg,ex3d(l,imin:imax,jmin:jmax,k),idx=idx,idy=idy)

            ex_z(ii,jj,kk)   = deriv1(zseg,ex3d(l,i,j,kmin:kmax),idx=idz)
            ex_xz(ii,jj,kk)  = deriv2(xseg,zseg,ex3d(l,imin:imax,j,kmin:kmax),idx=idx,idy=idz)
            ex_yz(ii,jj,kk)  = deriv2(yseg,zseg,ex3d(l,i,jmin:jmax,kmin:kmax),idx=idy,idy=idz)
            ex_xyz(ii,jj,kk) = deriv3(xseg,yseg,zseg,ex3d(l,imin:imax,jmin:jmax,kmin:kmax),idx=idx,idy=idy,idz=idz)

            ! Derivatives of the Y-component
            ey_x(ii,jj,kk)  = deriv1(xseg,ey3d(l,imin:imax,j,k),idx=idx)
            ey_y(ii,jj,kk)  = deriv1(yseg,ey3d(l,i,jmin:jmax,k),idx=idy)
            ey_xy(ii,jj,kk) = deriv2(xseg,yseg,ey3d(l,imin:imax,jmin:jmax,k),idx=idx,idy=idy)

            ey_z(ii,jj,kk)   = deriv1(zseg,ey3d(l,i,j,kmin:kmax),idx=idz)
            ey_xz(ii,jj,kk)  = deriv2(xseg,zseg,ey3d(l,imin:imax,j,kmin:kmax),idx=idx,idy=idz)
            ey_yz(ii,jj,kk)  = deriv2(yseg,zseg,ey3d(l,i,jmin:jmax,kmin:kmax),idx=idy,idy=idz)
            ey_xyz(ii,jj,kk) = deriv3(xseg,yseg,zseg,ey3d(l,imin:imax,jmin:jmax,kmin:kmax),idx=idx,idy=idy,idz=idz)

            ! Derivatives of the Z-component
            ez_x(ii,jj,kk)  = deriv1(xseg,ez3d(l,imin:imax,j,k),idx=idx)
            ez_y(ii,jj,kk)  = deriv1(yseg,ez3d(l,i,jmin:jmax,k),idx=idy)
            ez_xy(ii,jj,kk) = deriv2(xseg,yseg,ez3d(l,imin:imax,jmin:jmax,k),idx=idx,idy=idy)

            ez_z(ii,jj,kk)   = deriv1(zseg,ez3d(l,i,j,kmin:kmax),idx=idz)
            ez_xz(ii,jj,kk)  = deriv2(xseg,zseg,ez3d(l,imin:imax,j,kmin:kmax),idx=idx,idy=idz)
            ez_yz(ii,jj,kk)  = deriv2(yseg,zseg,ez3d(l,i,jmin:jmax,kmin:kmax),idx=idy,idy=idz)
            ez_xyz(ii,jj,kk) = deriv3(xseg,yseg,zseg,ez3d(l,imin:imax,jmin:jmax,kmin:kmax),idx=idx,idy=idy,idz=idz)

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
    evec = evec-dot_product(bhat,evec)*bhat
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
