import numpy as np
from scipy import linalg


def dipole_field(xv, b0=-31100):
    """
    """
    bx = lambda x,y,z: 3*b0*x*z/np.sqrt(x*x+y*y+z*z)**5
    by = lambda x,y,z: 3*b0*y*z/np.sqrt(x*x+y*y+z*z)**5
    bz = lambda x,y,z: b0*(2*z*z-x*x-y*y)/np.sqrt(x*x+y*y+z*z)**5
    
    bvec = np.r_[bx(*xv), by(*xv), bz(*xv)]
    
    return bvec


def dipole_gradient(xv, b0=-31100, r_e=6371, dx=0.01):
    """Estimate gradient of dipole field
    
    Uses finite differences becauseit saves dev time over
    calculating derivatives analytically.
    """
    
    bx = lambda x, y, z: 3*b0*x*z/np.sqrt(x*x+y*y+z*z)**5
    by = lambda x, y, z: 3*b0*y*z/np.sqrt(x*x+y*y+z*z)**5
    bz = lambda x, y, z: b0*(2*z*z-x*x-y*y)/np.sqrt(x*x+y*y+z*z)**5
    df = lambda f, x, d: (f(x+d)-f(x-d))/(2*d)
    
    dbxdx = df(lambda x: bx(x, xv[1], xv[2]), xv[0], dx)
    dbxdy = df(lambda x: bx(xv[0], x, xv[2]), xv[1], dx)
    dbxdz = df(lambda x: bx(xv[0], xv[1], x), xv[2], dx)
    gradbx = np.r_[dbxdx, dbxdy, dbxdz]/r_e
    
    dbydx = df(lambda x: by(x, xv[1], xv[2]), xv[0], dx)
    dbydy = df(lambda x: by(xv[0], x, xv[2]), xv[1], dx)
    dbydz = df(lambda x: by(xv[0], xv[1], x), xv[2], dx)
    gradby = np.r_[dbydx, dbydy, dbydz]/r_e
    
    dbzdx = df(lambda x: bz(x, xv[1], xv[2]), xv[0], dx)
    dbzdy = df(lambda x: bz(xv[0], x, xv[2]), xv[1], dx)
    dbzdz = df(lambda x: bz(xv[0], xv[1], x), xv[2], dx)
    gradbz = np.r_[dbzdx, dbzdy, dbzdz]/r_e
    
    gradb = np.c_[gradbx, gradby, gradbz]
    
    return gradb


def grad_bmag(xv, r_e=6371, dx=0.01):
    """Calculate the gradient of the magnetic field magnitude"""
    btot = lambda x, y, z: linalg.norm(dipole_field([x, y, z]))
    df = lambda f, x: (f(x + dx) - f(x - dx))/(2*dx)
    bgradx = df(lambda x: btot(x, xv[1], xv[2]), xv[0])
    bgrady = df(lambda x: btot(xv[0], x, xv[2]), xv[1])
    bgradz = df(lambda x: btot(xv[0], xv[1], x), xv[2])
    return np.r_[bgradx, bgrady, bgradz]/r_e
