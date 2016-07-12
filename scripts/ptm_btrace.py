"""
PTM_BTRACE

Magnetic field tracing and analysis for magnetic field data on scattered grids. This module is intended
to supercede the PTM_FIELDS_TRACING module.

This module defines a "GRIDDED_MAGNETIC_FIELD" object that reads in gridded magnetic field data
and allows for the fields to be arbitrarily evaluated, allows for field lines to be traced, and
allows for the determination of the magnetic equator.

GROUND_MAGNETIC_FIELD

Public Methods (see individual methods for documentation)
--configure_reader
--get_bhat
--get_bvec
--get_bmag
--find_min_B
--find_field_line
--trace_field_line
--trace_magnetic_equator

EXAMPLE USAGE

bfield=gridded_magnetic_field(istep=240,searchdir='/Users/jwoodroffe/Desktop/substorm_gridded/')
mageq=bfield.trace_magnetic_equator(6.6)

Jesse Woodroffe
6/28/2016
"""

from numpy import fromfile,pi,cos,sin,zeros,zeros_like,r_,c_,abs,linspace,size,argmin,dot,arccos,arctan2
from numpy.linalg import norm
from scipy.optimize import fsolve
from scipy.integrate import ode
from scipy.interpolate import RegularGridInterpolator

dtor = pi/180.0
rtod = 180.0/pi

class gridded_magnetic_field(object):
  """
  Base object for SWMF field line tracing
  """

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

    self.__get_grid()
    self.__get_fields()

  def __get_grid(self):
    """
    Read in spatial grids and set parameters
    """
    
    self.__xgrid=fromfile(self.__searchdir+'xgrid.bin')
    self.__ygrid=fromfile(self.__searchdir+'ygrid.bin')
    self.__zgrid=fromfile(self.__searchdir+'zgrid.bin')
    self.__nx=self.__xgrid.size
    self.__ny=self.__ygrid.size
    self.__nz=self.__zgrid.size

  def __get_fields(self):
    """
    Read in magnetic fields and create some inquiry functions.
    """
    
    bdims=(self.__nx,self.__ny,self.__nz)
    self.__bx=fromfile(self.__searchdir+'bx3d_%4.4i.bin' % self.__istep).reshape(bdims)
    self.__by=fromfile(self.__searchdir+'by3d_%4.4i.bin' % self.__istep).reshape(bdims)
    self.__bz=fromfile(self.__searchdir+'bz3d_%4.4i.bin' % self.__istep).reshape(bdims)
    bxi=RegularGridInterpolator((self.__xgrid,self.__ygrid,self.__zgrid),self.__bx)
    byi=RegularGridInterpolator((self.__xgrid,self.__ygrid,self.__zgrid),self.__by)
    bzi=RegularGridInterpolator((self.__xgrid,self.__ygrid,self.__zgrid),self.__bz)
    
    def get_bhat(s,xv):
      bhat=r_[bxi(xv),byi(xv),bzi(xv)]
      bhat/=norm(bhat)
      return bhat
        
    def get_bvec(xv):
      bvec=r_[bxi(xv),byi(xv),bzi(xv)]
      return bvec
        
    def get_bmag(xv):
      res=norm(r_[bxi(xv),byi(xv),bzi(xv)])
      return res   

    self.get_bhat = get_bhat
    self.get_bvec = get_bvec
    self.get_bmag = get_bmag

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
    
  def __find_tilt_axis(self,mlt=0.0,rad=4.0,lat=0.0):
    """
    Determine the local dipole tilt angles. That is, determine what the azimuthal and
    latitudinal angles are for the dipole magnetic moment based on the local fields. This is
    less accurate the further away you get from the dipole field. This routine is based on a
    simple vector algebra analysis of the general equation for a magnetic dipole.
    
    Input:
      mlt     float     optional    Magnetic local time in hours where evaluation should occur
      rad     float     optional    Radial distance in Earth radii where evaluation should occur
      lat     float     optional    Magnetic latitude in degrees where evaluation should occur
      
    Output
      tilt    float   Dipole tilt angle
      azim    float   Dipole azimuth angle
      mmag    float   Magnitude of the dipole magnetic moment
      
    If the field is a pure dipole, then the results of this routine should be independent of
    position. However, for a distorted field as is found in SWMF, this is probably not the case.
    
    Jesse Woodroffe
    6/27/2016
    """
    
    phi=self.__mlt_to_phi(mlt)
    x0=rad*cos(dtor*phi)*cos(dtor*lat)
    y0=rad*sin(dtor*phi)*cos(dtor*lat)
    z0=rad*sin(dtor*lat)
    pos=r_[x0,y0,z0]
    bvec=self.get_bvec(pos)
    r=norm(pos)
    
    rhat=pos/norm(pos)
    brad=dot(rhat,bvec)
    bhat=bvec/norm(bvec)
    bmag=norm(bvec)
    bsqr=bmag*bmag
    
    mpara=(1.5*brad**2-bsqr)/bmag
    mperp=1.5*brad*(bsqr*rhat-brad*bvec)/bsqr
    mtot=mperp+mpara*bhat
    mmag=norm(mtot)
    mhat=mtot/mmag

    tilt = rtod*arccos(mhat[2])
    azim = rtod*arctan2(mhat[1],mhat[0])
    mlt = self.__phi_to_mlt(azim)

    return tilt,azim,mmag*r**3
    
  def configure_reader(self,istep=None,searchdir=None):
    """
    Change data source
    """
    
    self.__init__(istep,searchdir)
   
  def find_min_B(self,mlt,rad,lat,ds=0.05):
    """
    Locate the position of minimum magnetic field along a given
    magnetic field line.
    
    Jesse Woodroffe
    6/28/2016
    """
    xt,yt,zt,bt,stot=self.trace_field_line(mlt,rad,lat,ds=ds)
    imin=argmin(bt)
    pos=r_[xt[imin],yt[imin],zt[imin]]
    return pos
    
  def find_field_line(self,mlt,r0,lat,ds=0.05):
    """
    Find the magnetic field line that has minimum B at a given radial distance and
    magnetic local time.
    
    Jesse Woodroffe
    6/28/2016
    """
    
    def rtfun(x):
      pos=self.find_min_B(mlt,r0,x,ds=ds)
      rval=norm(pos)
      sol=rval-r0
      return sol

    minlat=fsolve(rtfun,lat,xtol=1e-3)
    
    pos = self.find_min_B(mlt,r0,minlat,ds=ds)
    
    return pos

  def trace_field_line(self,mlt,rad,lat,ds=0.05,istep_max=100000):
    """
    Given a point in MLT, radial distance, latitude space, trace the corresponding
    magnetic field line to both ionospheres. Also determines the total length of the
    field line which may be useful e.g. in determining resonant frequencies.
  
    Jesse Woodroffe 5/18/16
    """
  
    phi=self.__mlt_to_phi(mlt)
    x0=rad*cos(dtor*phi)*cos(dtor*lat)
    y0=rad*sin(dtor*phi)*cos(dtor*lat)
    z0=rad*sin(dtor*lat)
  
    xp=zeros([istep_max])
    yp=zeros_like(xp)
    zp=zeros_like(xp)
    bp=zeros_like(xp)
    xm=zeros_like(xp)
    ym=zeros_like(xp)
    zm=zeros_like(xp)
    bm=zeros_like(xp)
  
    r=ode(self.get_bhat).set_integrator('vode',method='adams')
  
    stot = 0.0
  
    # Trace north
  
    ics=r_[x0,y0,z0]
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
            
        if(norm(r.y)<=1.0):
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
  
    stot+=abs(r.t)
  
    # Trace south
  
    ics=r_[x0,y0,z0]
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
      
        if(norm(r.y)<=1.0):
            doIntegrate=False
          
        if(r.y[0]<self.__xgrid.min()+1.0 or r.y[0]>self.__xgrid.max()-1.0):
            doIntegrate=False
          
        if(r.y[1]<self.__ygrid.min()+1.0 or r.y[1]>self.__ygrid.max()-1.0):
            doIntegrate=False
          
        if(r.y[2]<self.__zgrid.min()+1.0 or r.y[2]>self.__zgrid.max()-1.0):
            doIntegrate=False
  
    stot+=abs(r.t)
  
    xl = xm[:istep]
    yl = ym[:istep]
    zl = zm[:istep]
    bl = bm[:istep]
  
    xr=r_[xl[::-1][:-1],xk]
    yr=r_[yl[::-1][:-1],yk]
    zr=r_[zl[::-1][:-1],zk]
    br=r_[bl[::-1][:-1],bk]

    return xr,yr,zr,br,stot
   
  def trace_magnetic_equator(self,r0,ltmin=0.0,ltmax=23.0,nlt=24,ds=0.05):
    """
    Find all points at a given radial distance in the magnetic equator for a selected set of
    magnetic local times.

    Parallelization of this routine could produce speedup, but there are inherent difficulties in using
    the multiprocessing library (it won't pickle methods)
    
    Jesse Woodroffe
    6/28/2016
    """
    
    ltvec=linspace(ltmin,ltmax,nlt)
    xeq=zeros([nlt])
    yeq=zeros_like(xeq)
    zeq=zeros_like(xeq)
    dipoleTilt,azim,mag=self.__find_tilt_axis()
    for i,lt in enumerate(ltvec):
      phi=self.__mlt_to_phi(lt)
      if(dipoleTilt<90.0):
        lat0=dipoleTilt*cos(dtor*(phi-180.0))
      else:
        lat0=(180.0-dipoleTilt)*cos(dtor*(phi-180.0))
      pos = self.find_field_line(lt,r0,lat0,ds=ds)
      xeq[i]=pos[0]
      yeq[i]=pos[1]
      zeq[i]=pos[2]
    
    return c_[xeq,yeq,zeq]