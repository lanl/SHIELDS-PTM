"""
Generic orbital specification for particle tracing.
"""

import numpy as np
from scipy import integrate

__dtor = np.pi/180.0
__rtod = 180.0/np.pi

def orbit_eq(th,r_perigee,r_apogee,orbital_phase=0.0,inDegrees=True,tilt_deg=0.0):
    """
    Basic orbital equation in GSM equatorial plane, optionally rotated into "SM" coordinates

    Inputs
        th              Angle of orbit relative to the line of apsides
        r_perigee       Perigee distance in re
        r_apogee        Apogee distance in re
        orbital_phase   Optional (default 0) Angle between line of apsides and midnight
        inDegrees       Optional (default T) Is th in degrees?
        tilt_deg        Optional (default 0) SM tilt angle in degrees

    Outputs
        x               X position (re) in GSM/SM coordinates
        y               Y position (re) in GSM/SM coordinates
        z               Z position (re) in GSM/SM coordinates

    Author
        Jesse Woodroffe
        9/27/2016
    """

    GM = 1.541397e-06
    RE = 6371.0
    rho = r_apogee/r_perigee
    e=(rho-1.0)/(rho+1.0)
    phase = th+orbital_phase

    if(inDegrees):
        phase*=__dtor

    r=r_perigee*(1+e)/(1+e*np.cos(phase+__dtor*orbital_phase))
    thdot=np.sqrt(r_apogee*(1-e)*GM/r**4)
    v=r*RE*thdot

    x=r*np.cos(phase)
    y=r*np.sin(phase)

    # Rotate coordinates into the X-Z plane
    xprime=x*np.cos(__dtor*tilt_deg)
    zprime=x*np.sin(__dtor*tilt_deg)

    return xprime,y,zprime,v

def make_ephemeris(r_perigee,r_apogee,orbital_phase=0.0,nsamp=100,t0=0.0,inDegrees=True,tilt_deg=0.0):
    """
    Create a synthetic ephemeris based on a specified orbit.
    """

    rho = r_apogee/r_perigee
    e=(rho-1.0)/(rho+1.0)
    GM = 1.541397e-06
    tcof = np.sqrt((r_apogee*(1-e))**3/GM)

    th=np.linspace(0.0,360.0,nsamp)
    x,y,z,v=orbit_eq(th,r_perigee,r_apogee,orbital_phase,inDegrees,tilt_deg)

    # Given tilt_deg and position, we can determine the local direction of the velocity
    t = np.zeros_like(v)
    vx = np.zeros_like(v)
    vy = np.zeros_like(v)
    vz = np.zeros_like(v)

    zhat = np.r_[np.sin(__dtor*tilt_deg),0.0,np.cos(__dtor*tilt_deg)]

    for i in xrange(v.size):
        rhat=np.r_[x[i],y[i],z[i]]/np.sqrt(x[i]**2+y[i]**2+z[i]**2)
        vhat=np.cross(rhat,zhat)
        vx[i] = v[i]*vhat[0]
        vy[i] = v[i]*vhat[1]
        vz[i] = v[i]*vhat[2]
        t[i] = tcof*integrate.quad(lambda x:1.0/(1.0+e*np.cos(x+__dtor*orbital_phase))**2,__dtor*th[0],__dtor*th[i])[0]

    return np.c_[t,x,y,z,vx,vy,vz]

def write_ephemeris(filename,ephemeris):

    with open(filename,'w') as f:
        # Write a simple one-line header
        f.write('#{:11} {:12} {:12} {:12}'.format(*['Time(s)','X(Re)','Y(Re)','Z(Re)\n']))
        for i in xrange(ephemeris.shape[0]):
            dataline = '{:12.5f} {:12.7f} {:12.7f} {:12.7f}\n'.format(*ephemeris[i,:4])
            f.write(dataline)