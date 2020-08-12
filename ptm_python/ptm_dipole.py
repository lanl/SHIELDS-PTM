import numpy as np
from scipy import linalg

def dipole_field(xv):
    b0 = -31100.0
    bx = lambda x,y,z:3*b0*x*z/np.sqrt(x*x+y*y+z*z)**5
    by = lambda x,y,z:3*b0*y*z/np.sqrt(x*x+y*y+z*z)**5
    bz = lambda x,y,z:b0*(2*z*z-x*x-y*y)/np.sqrt(x*x+y*y+z*z)**5
    
    bvec = np.r_[bx(*xv),by(*xv),bz(*xv)]
    
    return bvec

def dipole_gradient(xv):
    # Estimate gradient of dipole field using finite differences
    # because, honestly, I don't feel like doing all the derivatives
    # analytically.
    
    b0 = -31100.0
    bx = lambda x,y,z:3*b0*x*z/np.sqrt(x*x+y*y+z*z)**5
    by = lambda x,y,z:3*b0*y*z/np.sqrt(x*x+y*y+z*z)**5
    bz = lambda x,y,z:b0*(2*z*z-x*x-y*y)/np.sqrt(x*x+y*y+z*z)**5
    df = lambda f,x,d:(f(x+d)-f(x-d))/(2*d)
    
    dbxdx = df(lambda x:bx(x,xv[1],xv[2]),xv[0],0.01)
    dbxdy = df(lambda x:bx(xv[0],x,xv[2]),xv[1],0.01)
    dbxdz = df(lambda x:bx(xv[0],xv[1],x),xv[2],0.01)
    gradbx = np.r_[dbxdx,dbxdy,dbxdz]/6371
    
    dbydx = df(lambda x:by(x,xv[1],xv[2]),xv[0],0.01)
    dbydy = df(lambda x:by(xv[0],x,xv[2]),xv[1],0.01)
    dbydz = df(lambda x:by(xv[0],xv[1],x),xv[2],0.01)
    gradby = np.r_[dbydx,dbydy,dbydz]/6371
    
    dbzdx = df(lambda x:bz(x,xv[1],xv[2]),xv[0],0.01)
    dbzdy = df(lambda x:bz(xv[0],x,xv[2]),xv[1],0.01)
    dbzdz = df(lambda x:bz(xv[0],xv[1],x),xv[2],0.01)
    gradbz = np.r_[dbzdx,dbzdy,dbzdz]/6371
    
    gradb=np.c_[gradbx,gradby,gradbz]
    
    return gradb

def grad_bmag(xv):
    # Calculate the gradient of the magnetic field magnitude
    d=0.01
    btot = lambda x,y,z:linalg.norm(dipole_field([x,y,z]))
    df = lambda f,x:(f(x+d)-f(x-d))/(2*d)
    bgradx = df(lambda x:btot(x,xv[1],xv[2]),xv[0])
    bgrady = df(lambda x:btot(xv[0],x,xv[2]),xv[1])
    bgradz = df(lambda x:btot(xv[0],xv[1],x),xv[2])
    return np.r_[bgradx,bgrady,bgradz]/6371