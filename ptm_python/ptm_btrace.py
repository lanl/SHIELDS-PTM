"""
PTM_BTRACE

Magnetic field tracing and analysis for magnetic field data on scattered grids. This module is intended
to supercede the PTM_FIELDS_TRACING module.

This module defines a "GRIDDED_MAGNETIC_FIELD" object that reads in gridded magnetic field data
and allows for the fields to be arbitrarily evaluated, allows for field lines to be traced, and
allows for the determination of the magnetic equator.

GRIDDED_MAGNETIC_FIELD

Public Methods (see individual methods for documentation)

configure_reader
get_bhat
get_bvec
get_bmag
find_min_B
find_field_line
trace_field_line
trace_magnetic_equator

EXAMPLE USAGE

import ptm_btrace as bt
bfield=bt.gridded_magnetic_field(istep=360,searchdir='/Users/jwoodroffe/Workspace/Projects/PTM/Events/7-18-2013/gridded/')
mageq=bfield.trace_magnetic_equator(6.6)

Jesse Woodroffe
6/28/2016
"""

import numpy as np
from numpy import linalg
from scipy import optimize
from scipy import integrate
from scipy import interpolate

class gridded_magnetic_field(object):
  """
  Base object for SWMF field line tracing
  """

  __dtor = np.pi/180.0
  __rtod = 180.0/np.pi

  def __init__(self,istep=None,searchdir=None):

    if(istep==None):
      self.__istep=0
    else:
      self.__istep=istep

    if(any([searchdir==None,searchdir=='',searchdir=='.'])):
      self.__searchdir=''
    else:
      if(searchdir[-1]!='/'):
        self.__searchdir=searchdir+'/'
      else:
        self.__searchdir=searchdir

    self.__xgrid=np.fromfile(self.__searchdir+'xgrid.bin')
    self.__ygrid=np.fromfile(self.__searchdir+'ygrid.bin')
    self.__zgrid=np.fromfile(self.__searchdir+'zgrid.bin')
    self.__nx=self.__xgrid.size
    self.__ny=self.__ygrid.size
    self.__nz=self.__zgrid.size

    bdims=(self.__nx,self.__ny,self.__nz)

    self.__bx=np.fromfile(self.__searchdir+'bx3d_%4.4i.bin' % self.__istep).reshape(bdims)
    self.__by=np.fromfile(self.__searchdir+'by3d_%4.4i.bin' % self.__istep).reshape(bdims)
    self.__bz=np.fromfile(self.__searchdir+'bz3d_%4.4i.bin' % self.__istep).reshape(bdims)

    bxi=interpolate.RegularGridInterpolator((self.__xgrid,self.__ygrid,self.__zgrid),self.__bx)
    byi=interpolate.RegularGridInterpolator((self.__xgrid,self.__ygrid,self.__zgrid),self.__by)
    bzi=interpolate.RegularGridInterpolator((self.__xgrid,self.__ygrid,self.__zgrid),self.__bz)

    def bhat_ode(s,xv):
      bhat=np.r_[bxi(xv),byi(xv),bzi(xv)]
      bhat/=linalg.norm(bhat)
      return bhat

    def bhat_odeint(xv,s):
      bhat=np.r_[bxi(xv),byi(xv),bzi(xv)]
      bhat/=linalg.norm(bhat)
      return bhat

    def get_bhat(xv):
      bhat=np.r_[bxi(xv),byi(xv),bzi(xv)]
      bhat/=linalg.norm(bvec)
      return bhat

    def get_bvec(xv):
      bvec=np.r_[bxi(xv),byi(xv),bzi(xv)]
      return bvec

    def get_bmag(xv):
      res=linalg.norm(np.r_[bxi(xv),byi(xv),bzi(xv)])
      return res

    # Create persistent methods
    self.get_bhat = get_bhat
    self.get_bvec = get_bvec
    self.get_bmag = get_bmag
    self.__bhat_ode = bhat_ode
    self.__bhat_odeint = bhat_odeint

  def __mlt_to_phi(self,myMlt):
    # Convert from magnetic local time in hours to azimuthal angle in degrees
    if(myMlt < 12.0):
      mlt=myMlt+24
    else:
      mlt=myMlt

    res = 15*(mlt-12.0)

    return res

  def __phi_to_mlt(self,myPhi):
    # Convert from azimuthal angle in degrees to magnetic local time in hours
    mlt = myPhi/15.0+12
    if(mlt >= 24.0): mlt-=24

    return mlt

  def configure_reader(self,istep=None,searchdir=None):
    """
    Change data source
    """

    self.__init__(istep,searchdir)

  def trace_field_line_section(self,mlt,rad,lat,ds=0.05,smax=2.0):
    """
    Given a point in MLT, radial distance, latitude space, trace a section of the
    corresponding field line to a distance +-smax.

    Jesse Woodroffe
    7/27/2016
    """

    phi=self.__mlt_to_phi(mlt)
    x0=rad*np.cos(self.__dtor*phi)*cos(self.__dtor*lat)
    y0=rad*np.sin(self.__dtor*phi)*cos(self.__dtor*lat)
    z0=rad*np.sin(self.__dtor*lat)

    swant=np.arange(0.0,smax,ds)
    ics=r_[x0,y0,z0]

    yplus=integrate.odeint(self.__bhat_odeint,ics,swant)
    yminus=integrate.odeint(self.__bhat_odeint,ics,-swant)

    y=np.vstack([yminus[::-1,:],yplus[1:,:]])

    bt = np.array([self.get_bmag(y[i,:]) for i in enumerate(y[:,0])])

    return y[:,0], y[:,1], y[:,2], bt

  def trace_field_line(self,mlt,rad,lat,ds=0.05,istep_max=100000):
    """
    Given a point in MLT, radial distance, latitude space, trace the corresponding
    magnetic field line to both ionospheres. Also determines the total length of the
    field line which may be useful e.g. in determining resonant frequencies.

    Jesse Woodroffe
    5/18/16
    """

    phi=self.__mlt_to_phi(mlt)
    x0=rad*np.cos(self.__dtor*phi)*np.cos(self.__dtor*lat)
    y0=rad*np.sin(self.__dtor*phi)*np.cos(self.__dtor*lat)
    z0=rad*np.sin(self.__dtor*lat)

    xp=np.zeros([istep_max])
    yp=np.zeros_like(xp)
    zp=np.zeros_like(xp)
    bp=np.zeros_like(xp)
    xm=np.zeros_like(xp)
    ym=np.zeros_like(xp)
    zm=np.zeros_like(xp)
    bm=np.zeros_like(xp)

    r=integrate.ode(self.__bhat_ode).set_integrator('vode',method='adams')

    stot = 0.0

    # Trace north

    ics=np.r_[x0,y0,z0]
    s0=0.0

    xp[0]=x0
    yp[0]=y0
    zp[0]=z0
    bp[0]=self.get_bmag(ics)

    r.set_initial_value(ics,s0)
    doIntegrate = True
    istep=0

    while r.successful() and doIntegrate:
        r.integrate(r.t+ds)

        istep+=1

        xp[istep]=r.y[0]
        yp[istep]=r.y[1]
        zp[istep]=r.y[2]
        bp[istep]=self.get_bmag(r.y)

        if(linalg.norm(r.y)<=1.0):
            doIntegrate=False

        if(r.y[0]<self.__xgrid.min()+1.0 or r.y[0]>self.__xgrid.max()-1.0):
            doIntegrate=False

        if(r.y[1]<self.__ygrid.min()+1.0 or r.y[1]>self.__ygrid.max()-1.0):
            doIntegrate=False

        if(r.y[2]<self.__zgrid.min()+1.0 or r.y[2]>self.__zgrid.max()-1.0):
            doIntegrate=False

    xk=xp[:istep]
    yk=yp[:istep]
    zk=zp[:istep]
    bk=bp[:istep]

    stot+=np.abs(r.t)

    # Trace south

    ics=np.r_[x0,y0,z0]
    s0=0.0

    xm[0]=x0
    ym[0]=y0
    zm[0]=z0
    bm[0]=self.get_bmag(ics)

    r.set_initial_value(ics,s0)
    doIntegrate = True
    istep=0

    while r.successful() and doIntegrate:
        r.integrate(r.t-ds)

        istep+=1

        xm[istep]=r.y[0]
        ym[istep]=r.y[1]
        zm[istep]=r.y[2]
        bm[istep]=self.get_bmag(r.y)

        if(linalg.norm(r.y)<=1.0):
            doIntegrate=False

        if(r.y[0]<self.__xgrid.min()+1.0 or r.y[0]>self.__xgrid.max()-1.0):
            doIntegrate=False

        if(r.y[1]<self.__ygrid.min()+1.0 or r.y[1]>self.__ygrid.max()-1.0):
            doIntegrate=False

        if(r.y[2]<self.__zgrid.min()+1.0 or r.y[2]>self.__zgrid.max()-1.0):
            doIntegrate=False

    stot+=np.abs(r.t)

    xl = xm[:istep]
    yl = ym[:istep]
    zl = zm[:istep]
    bl = bm[:istep]

    xr=np.r_[xl[::-1][:-1],xk]
    yr=np.r_[yl[::-1][:-1],yk]
    zr=np.r_[zl[::-1][:-1],zk]
    br=np.r_[bl[::-1][:-1],bk]

    return xr,yr,zr,br,stot

  def find_min_B(self,mlt,rad,lat,ds=0.05):
    """
    Locate the position of minimum magnetic field along a given
    magnetic field line.

    Jesse Woodroffe
    6/28/2016
    """
    xt,yt,zt,bt,stot=self.trace_field_line(mlt,rad,lat,ds=ds)
    imin=np.argmin(bt)
    pos=np.r_[xt[imin],yt[imin],zt[imin]]
    return pos

  def __find_min_B_root(self,mlt,rad,lat,ds=0.05,smax=1.0):

    xt,yt,zt,bt=self.trace_field_line_section(mlt,rad,lat,ds=ds)
    imin=np.argmin(bt)
    pos=np.r_[xt[imin],yt[imin],zt[imin]]
    return pos

  def __find_field_line_root(self,mlt,r0,lat,ds=0.05,smax=1.0):
    """
    Find the magnetic field line that has minimum B at a given radial distance and
    magnetic local time.

    Jesse Woodroffe
    6/28/2016
    """

    def rtfun(x):
      pos=self.__find_min_B_root(mlt,r0,x,ds=ds,smax=smax)
      rval=linalg.norm(pos)
      sol=rval-r0
      return sol

    minlat=optimize.newton(rtfun,lat,tol=1e-3)

    phi = self.__mlt_to_phi(mlt)
    x = r0*np.cos(self.__dtor*minlat)*np.cos(self.__dtor*phi)
    y = r0*np.cos(self.__dtor*minlat)*np.sin(self.__dtor*phi)
    z = r0*np.sin(self.__dtor*minlat)

    pos = np.r_[x,y,z]
    return pos

  def __find_field_line_arc(self,r0,mlt,arcSize=15.0,numpts=100,rlev=2,smplane=False):
    """
    Find the location of minimum B along an arc at a given radial distance and magnetic local time.
    For simple magnetic geometries, this locates the minimum-B field line. For more complicated sitatuions,
    this provides a good guess for the more exhaustive __find_field_line_root. The combination of these
    two functionalities is provided by the general find_field_line routine below.

    Jesse Woodroffe
    7/27/2016
    """

    # Get initial guess from dipole approximation
    lat0,phi=self.__dipole_eq(mlt,r0)

    if not smplane: # Improve by searching along arc

      # Intial bounds for arc search
      lmin=lat0-arcSize
      lmax=lat0+arcSize

      for ilev in xrange(rlev):

        # Parameterize the arc
        latv=self.__dtor*np.linspace(lmin,lmax,numpts)

        xv=r0*np.cos(latv)*np.cos(self.__dtor*phi)
        yv=r0*np.cos(latv)*np.sin(self.__dtor*phi)
        zv=r0*np.sin(latv)

        bmat=np.zeros([latv.size,4])

        # Calculate magnetic field magnitudes along the arc
        for i in xrange(latv.size):
          bmat[i,:3]=self.get_bvec(r_[xv[i],yv[i],zv[i]])
          bmat[i, 3]=linalg.norm(bmat[i,:3])

        # Refine bracketing
        imin=np.argmin(bmat[:,-1])
        isub=imin-2 if imin>2 else 0
        iadd=imin+2 if imin<numpts-3 else numpts-1
        lmin=self.__rtod*latv[isub]
        lmax=self.__rtod*latv[iadd]

      xr,yr,zr = xv[imin], yv[imin], zv[imin]

    else:

      xr=r0*np.cos(self.__dtor*lat0)*np.cos(self.__dtor*phi)
      yr=r0*np.cos(self.__dtor*lat0)*np.sin(self.__dtor*phi)
      zr=r0*np.sin(self.__dtor*lat0)

    return xr, yr, zr

  def __dipole_eq(self,mlt,r):
    """
    Determine the location of the specified point on the dipole magnetic equator.

    Jesse Woodroffe
    7/28/2016
    """

    # Calculate the magnetic moment vector
    bvec=self.get_bvec(np.r_[3.0,0.0,0.0])
    rhat=np.r_[1.0,0.0,0.0]
    mvec=1.5*dot(rhat,bvec)*rhat-bvec

    # Determine the planar coefficients for the magnetic equator
    mhat=mvec/linalg.norm(mvec)
    a=mhat[0]
    b=mhat[1]
    c=mhat[2]

    # Determine the location by enforcing planar and fixed distance constraints
    phi = self.__mlt_to_phi(mlt)*self.__dtor
    rho=r*np.abs(c)/np.sqrt((a*np.cos(phi)+b*np.sin(phi))**2+c**2)
    x=rho*np.cos(phi)
    y=rho*np.sin(phi)
    z=r*(a*np.cos(phi)+b*np.sin(phi))/sqrt((a*np.cos(phi)+b*np.sin(phi))**2+c**2)

    # Return the data
    lat=np.arcsin(z/r)*self.__rtod
    azm=phi*self.__rtod

    return lat,azm

  def find_field_line(self,r0,mlt,arcSize=15.0,numpts=100.0,rlev=2,smplane=False,refine=False,refine_ds=0.05,smax=1.0):
    """
    Locate a specified field line using an arc-based B-minimizer and rootfinding-based optimizer to refine.

    Jesse Woodroffe
    7/27/2016
    """
    pos=self.__find_field_line_arc(r0,mlt,arcSize=arcSize,numpts=numpts,rlev=rlev,smplane=smplane)
    if(refine):
      lat=self.__rtod*np.arcsin(pos[2]/linalg.norm(pos))
      try:
        pos=self.__find_field_line_root(mlt,r0,lat,ds=refine_ds,smax=smax)
      except: # Occasionally higher resolution is required to get a convergent solution
        pos=self.__find_field_line_root(mlt,r0,lat,ds=0.01,smax=smax)
    return pos

  def trace_magnetic_equator(self,r0,ltmin=0.0,ltmax=23.0,nlt=24,ds=0.05,refine=False,smplane=False):
    """
    Find all points at a given radial distance in the magnetic equator for a selected set of
    magnetic local times.

    Jesse Woodroffe
    6/28/2016
    """

    ltvec=np.linspace(ltmin,ltmax,nlt)
    xeq=np.zeros([nlt])
    yeq=np.zeros_like(xeq)
    zeq=np.zeros_like(xeq)
    for i,lt in enumerate(ltvec):
      pos = self.find_field_line(r0,lt,refine_ds=ds,refine=refine,smplane=smplane)
      xeq[i]=pos[0]
      yeq[i]=pos[1]
      zeq[i]=pos[2]

    return np.c_[xeq,yeq,zeq]
