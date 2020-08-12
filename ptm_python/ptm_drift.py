import numpy as np
from scipy import optimize
from scipy import integrate
from scipy import interpolate

def wc(q,mc2,B):
    """
    Calculate the electron gyrofrequency
    
    q is the charge in multiples of the fundamental
    mc2 is the particle rest mass in MeV
    B is the magnetic field magnitude in nT
    """
    
    cof=89.8755311
    res=cof*q*B/mc2
    return res

def T(alpha):
    # Auxiliary integral used in drift calculations
    # Hamlin 1961 Equation 16
    mu=np.sin(np.radians(alpha))
    thetac=optimize.brentq(lambda x:np.sin(x)**6-mu*mu*np.sqrt(1+3*np.cos(x)**2),0,np.pi/2)
    num=lambda x:np.sin(x)*np.sqrt(1+3*np.cos(x)**2)
    den=lambda x:np.sqrt(1-mu*mu*np.sqrt(1+3*np.cos(x)**2)/np.sin(x)**6)
    res=integrate.quad(lambda x:num(x)/den(x),thetac,np.pi/2)[0]
    return res

def E(alpha):
    # Auxiliary integral used in drift calculations    
    # Hamlin 1961 Equation 22 with correction from corrigendum
    mu=np.sin(np.radians(alpha))
    thetac=optimize.brentq(lambda x:np.sin(x)**6-mu*mu*np.sqrt(1+3*np.cos(x)**2),0,np.pi/2)
    num=lambda x:(np.sin(x)**3*(1+np.cos(x)**2))*(1-0.5*mu*mu*np.sqrt(1+3*cos(x)**2)/np.sin(x)**6)
    den=lambda x:(1+3*np.cos(x)**2)**(3/2)*np.sqrt(1-mu*mu*np.sqrt(1+3*cos(x)**2)/np.sin(x)**6)
    res=integrate.quad(lambda x:num(x)/den(x),thetac,np.pi/2)[0]
    return res

def T_drift(Ekin,alpha,L,q=1,mc2=0.511):
    """
    Calculate the drift period of a particle in the Earth's dipole field
    Based on Hamlin 1961 Equation 21
    
    Ekin is the particle kinetic energy in MeV
    alpha is the pitch angle in degrees
    L is the dipole L-shell parameter
    q is the charge in multiples of the fundamental
    mc2 is the particle rest mass in MeV (e.g., 0.511 for electrons, 938 for protons)
    """
    
    b0=3.11e4/L**3
    r0=6371e3*L
    wc=10*(2.9979248)**2*q*b0/mc2
    if alpha > 89:
        ratio=0.35+0.15*np.sin(np.radians(alpha))
    else:
        ratio=E(alpha)/T(alpha)
    gam=1+Ekin/mc2
    v=2.9979248e8*np.sqrt(gam*gam-1)/gam
    wd=3*(gam/wc)*(v/r0)**2*ratio
    return 2*np.pi/wd