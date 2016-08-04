#!/usr/bin/python

"""
When run as a command-line script, traces the magnetic equator of a specified magnetic field model
that is described in the files b{x,y,z}3d_xxxx.bin. By default, the equator is calculated at 24 
points at geosynchronous distances by searching for the field line whose minimum  B occurs at that
MLT.

This file contains routines to trace magnetic field lines, determine the location of minimum B along a 
magnetic field line, to determine the equatorial location of a minimum B location for a given MLT and 
radial distance, and to map the equator at a given radial distance for a specified range of MLTs.

Example usage:
import ptm_fields_tracing
from pylab import plot,show
from numpy import pi,arcsin
from numpy.linalg import norm

ptm_fields_tracing.get_B_fields(60)
eqPos=ptm_fields_tracing.find_field_line(0.0,6.6,23.0,dlat=0.25,verbose=False)
eqLat=180/pi*arcsin(eqPos[2]/norm(eqPos))
xr,yr,zr=ptm_fields_tracing.trace_field_line(0.0,6.6,eqLat)
xm,ym,zm,pos=ptm_fields_tracing.find_min_B_position(0.0,6.6,eqLat)
plot(xr,zr)
plot(pos[0],pos[2],'ro',ms=8)
show()

(There are usually some warning messages about integration accuracy, but the results are typically good.)

Jesse Woodroffe
5/23/2016
"""

from numpy import fromfile,pi,cos,sin,zeros,zeros_like,r_,c_,abs,linspace,size
from numpy.linalg import norm
from scipy.integrate import ode
from scipy.interpolate import RegularGridInterpolator
from sys import argv

global nx, ny, nz
global xgrid, ygrid, zgrid
global bx, by, bz

def __get_grid(searchDir=''):
  """
  The searchDir argument specifies which directory to look in for grid files
  """

  global nx, ny, nz
  global xgrid, ygrid, zgrid

  xgrid=fromfile(searchDir+'xgrid.bin')
  ygrid=fromfile(searchDir+'ygrid.bin')
  zgrid=fromfile(searchDir+'zgrid.bin')

  nx = xgrid.size
  ny = ygrid.size
  nz = zgrid.size

  return
  
def get_B_fields(istep,searchDir=''):

    """
    The integer istep argument specifies which timestep to read in.
    The searchDir argument specifies which directory to look in for magnetic field and grid files
    """

    global nx, ny, nz
    global bx, by, bz

    # Handle a couple of special cases
    if(len(searchDir)>0):
      if(searchDir=='.'):
        # User specified pwd using dot; strip the dot 
        searchDir=''
      elif(searchDir[-1]!='/'): 
        # User forgot to add a slash at the end of the directory name; add a slash
        searchDir+='/'

    __get_grid(searchDir=searchDir)
      
    bx=fromfile(searchDir+'bx3d_%4.4i.bin' % istep).reshape((nx,ny,nz))
    by=fromfile(searchDir+'by3d_%4.4i.bin' % istep).reshape((nx,ny,nz))
    bz=fromfile(searchDir+'bz3d_%4.4i.bin' % istep).reshape((nx,ny,nz))

    return
    
def find_min_B_position(mlt,rad,lat,ds=0.01,istep_max=10000):
    """
    Given a point in MLT, radial distance, latitude space, determine the Cartesian
    position of the minimum-B point (aka the magnetic equator) for that field line
    
    Jesse Woodroffe 5/18/16
    """
   
    phi=pi*(1.0-mlt/12.0)
    x0=rad*cos(phi)*cos(pi*lat/180.)
    y0=rad*sin(phi)*cos(pi*lat/180.)
    z0=rad*sin(pi*lat/180.)
    
    bxi = RegularGridInterpolator((xgrid,ygrid,zgrid),bx)
    byi = RegularGridInterpolator((xgrid,ygrid,zgrid),by)
    bzi = RegularGridInterpolator((xgrid,ygrid,zgrid),bz)
    
    xtr=zeros([istep_max])
    ytr=zeros_like(xtr)
    ztr=zeros_like(xtr)
    
    def dbfun(s,xv):
        bhat=r_[bxi(xv),byi(xv),bzi(xv)]; bhat/=norm(bhat)
        return bhat
    
    def bmag(xv):
        res=norm(r_[bxi(xv),byi(xv),bzi(xv)])
        return res
    
    r=ode(dbfun).set_integrator('vode',method='adams')
    ics=r_[x0,y0,z0]
    bold=bmag(ics)
    s0=0.0
    r.set_initial_value(ics,s0)
    doIntegrate = True
    istep=0
    
    xtr[0]=x0
    ytr[0]=y0
    ztr[0]=z0
    
    while r.successful() and doIntegrate:
        yold=r.y.copy()
        r.integrate(r.t+ds)
        bnow=bmag(r.y)
        if(bnow>bold):
            if(istep==0): # Integration is proceeding away from equator
                ds*=-1
            else: # We've just passed the minimum, pass back the last point
                doIntegrate = False
        istep+=1
        if(istep>istep_max):
            doIntegrate=False
        xtr[istep]=r.y[0]
        ytr[istep]=r.y[1]
        ztr[istep]=r.y[2]
        bold=bnow.copy()

    return xtr[:istep],ytr[:istep],ztr[:istep],yold
    
def find_field_line(mlt,r0,lat0,dlat=5.0,verbose=True,epserr=1e-2):
    """
    Given a particular mlt and radial distance, search for the magnetic field line
    that actually has minimum B with those properties.
    
    Jesse Woodroffe
    5/18/2016
    """
    
    xtm,ytm,ztm,yo=find_min_B_position(mlt,r0,lat0)
    xtu,ytu,ztu,yu=find_min_B_position(mlt,r0,lat0+dlat)
   
    # Determine initial search direction
    if(norm(yo) < norm(yu)):
        dlat*=-1.0
        lat = lat0
        pos=yo
    else:
        pos=yu
        lat = lat0+dlat

    doSearch = True
    isearch = 0
        
    while doSearch:
        
        pold = pos.copy()

        lat+=dlat
        xt,yt,zt,pos=find_min_B_position(mlt,r0,lat)

        isearch+=1

        if(verbose):
            print lat, norm(pos), norm(pold)
        
        if(norm(pos) > norm(pold) or isearch>100 or abs(r0-norm(pos))<epserr):
            doSearch = False
            
    return pos

def trace_field_line(mlt,rad,lat,ds=0.01,istep_max=100000):
    """
    Given a point in MLT, radial distance, latitude space, trace the corresponding
    magnetic field line to both ionospheres. Also determines the total length of the
    field line which may be useful e.g. in determining resonant frequencies.
    
    Jesse Woodroffe 5/18/16
    """
    
    phi=pi*(1.0-mlt/12.0)
    x0=rad*cos(phi)*cos(pi*lat/180.)
    y0=rad*sin(phi)*cos(pi*lat/180.)
    z0=rad*sin(pi*lat/180.)
       
    bxi = RegularGridInterpolator((xgrid,ygrid,zgrid),bx)
    byi = RegularGridInterpolator((xgrid,ygrid,zgrid),by)
    bzi = RegularGridInterpolator((xgrid,ygrid,zgrid),bz)
    
    xp=zeros([istep_max])
    yp=zeros_like(xp)
    zp=zeros_like(xp)
    xm=zeros_like(xp)
    ym=zeros_like(xp)
    zm=zeros_like(xp)
    
    def dbfun(s,xv):
        bhat=r_[bxi(xv),byi(xv),bzi(xv)]; bhat/=norm(bhat)
        return bhat
    
    r=ode(dbfun).set_integrator('vode',method='adams')
    
    stot = 0.0
    
    # Trace north
    
    ics=r_[x0,y0,z0]
    s0=0.0
    
    r.set_initial_value(ics,s0)
    doIntegrate = True
    istep=0
    
    xp[0]=x0
    yp[0]=y0
    zp[0]=z0
        
    while r.successful() and doIntegrate:
        r.integrate(r.t+ds)

        istep+=1
        
        xp[istep]=r.y[0]
        yp[istep]=r.y[1]
        zp[istep]=r.y[2]
              
        if(norm(r.y)<=1.0):
            doIntegrate=False
            
        if(r.y[0]<xgrid.min()+1.0 or r.y[0]>xgrid.max()-1.0):
            doIntegrate=False
            
        if(r.y[1]<ygrid.min()+1.0 or r.y[1]>ygrid.max()-1.0):
            doIntegrate=False
            
        if(r.y[2]<zgrid.min()+1.0 or r.y[2]>zgrid.max()-1.0):
            doIntegrate=False
    
    xk=xp[:istep]
    yk=yp[:istep]
    zk=zp[:istep]
    
    stot+=abs(r.t)
    
    # Trace south
    
    ics=r_[x0,y0,z0]
    s0=0.0
    
    r.set_initial_value(ics,s0)
    doIntegrate = True
    istep=0
    
    xm[0]=x0
    ym[0]=y0
    zm[0]=z0
    
    while r.successful() and doIntegrate:
        r.integrate(r.t-ds)

        istep+=1
        
        xm[istep]=r.y[0]
        ym[istep]=r.y[1]
        zm[istep]=r.y[2]
        
        if(norm(r.y)<=1.0):
            doIntegrate=False
            
        if(r.y[0]<xgrid.min()+1.0 or r.y[0]>xgrid.max()-1.0):
            doIntegrate=False
            
        if(r.y[1]<ygrid.min()+1.0 or r.y[1]>ygrid.max()-1.0):
            doIntegrate=False
            
        if(r.y[2]<zgrid.min()+1.0 or r.y[2]>zgrid.max()-1.0):
            doIntegrate=False
    
    stot+=abs(r.t)
    
    xl = xm[:istep]
    yl = ym[:istep]
    zl = zm[:istep]
    
    xr=r_[xl[::-1][:-1],xk]
    yr=r_[yl[::-1][:-1],yk]
    zr=r_[zl[::-1][:-1],zk]
  
    return xr,yr,zr,stot
    
def trace_magnetic_equator(r0=6.6,ltmin=0.0,ltmax=23.0,nlt=24,dipoleTilt=20.0):
  ltvec=linspace(ltmin,ltmax,nlt)
  xeq=zeros([nlt])
  yeq=zeros_like(xeq)
  zeq=zeros_like(xeq)
  for i,lt in enumerate(ltvec):
    lat0=dipoleTilt*cos(pi*lt/12.0)-5*(1-cos(pi*lt/12.0))
    pos = find_field_line(lt,r0,lat0,dlat=1.0,verbose=False)
    xeq[i]=pos[0]
    yeq[i]=pos[1]
    zeq[i]=pos[2]
    
  return c_[xeq,yeq,zeq]
    
if __name__ == "__main__":
    index=argv[1]
    get_B_fields(index)
    mageq=trace_magnetic_equator()
    
    # Write equator data to file
    eqfile=open('mageq_%4.4i.dat' % index,'w')
    for i in xrange(size(mageq,0)):
      eqfile.writeline('%12.5f %12.5f %12.5f\n' % mageq[i,:])
    eqfile.close()