"""
This file contains an object-oriented implementation of the Tsyganenko-Mukai
plasma sheet model. This model is based off of a statistical analysis of
measurements from the GEOTAIL spacecraft. See description in source paper for
more information.

References:

Tsyganenko, N. A., and T. Mukai (2003), Tail plasma sheet models derived from Geotail
particle data, J. Geophys. Res., 108, 1136, doi:10.1029/2002JA009707, A3

"""

import numpy as np

class plasma_sheet(object):

    __aT=np.array([ 0.0000, 1.6780,-0.1606, 1.6690, 4.8200, 2.8550,-0.6020,-0.8360,
                 -2.4910, 0.2568, 0.2249, 0.1887,-0.4458,-0.0331,-0.0241,-2.6890,
                  1.2220])
    __aN=np.array([ 0.0000,-0.1590, 0.6080, 0.5055, 0.0796, 0.2746, 0.0361,-0.0342,
                 -0.7935, 1.1620, 0.4756, 0.7117])
    __aP=np.array([ 0.0000, 0.0570, 0.5240, 0.0908, 0.5270, 0.0780,-4.4220,-1.5330,
                 -1.2170, 2.5400, 0.3200, 0.7540, 1.0480,-0.0740, 1.0150])

    __bnorm = 5.0
    __vnorm = 500.0
    __nnorm = 10.0
    __rnorm = 10.0
    __pnorm = 3.0
    __dtor = np.pi/180.0
    __rtod = 180.0/np.pi

    def __init__(self,bperp=__bnorm,theta=90.0,vx=__vnorm,n=__nnorm,p=__pnorm):

        self.bperp = bperp/self.__bnorm
        self.theta = theta
        self.vsw = vx/self.__vnorm
        self.nsw = n/self.__nnorm
        self.psw = p/self.__pnorm
        self.bz = self.bperp*np.cos(self.__dtor*theta)

        if(self.bz > 0.0):
            self.bzs = 0.0
            self.bzn = self.bz
        else:
            self.bzn = 0.0
            self.bzs = -self.bz

        self.fsw = self.bperp*np.sqrt(np.sin(0.5*theta*self.__dtor))

        return

    def set_parameters(self,bperp=None,theta=None,vx=None,n=None,p=None):

        if not bperp is None:
            self.bperp = bperp/self.__bnorm

        if not theta is None:
            self.theta = theta

        if not vx is None:
            self.vsw = vx/self.__vnorm

        if not n is None:
            self.nsw = n/self.__nnorm

        if not p is None:
            self.psw = p/self.__pnorm

        if (not bperp is None) or (not theta is None):
            self.bz = self.bperp*np.cos(self.__dtor*self.theta)
            if self.bz > 0.0:
                self.bzs = 0.0
                self.bzn = self.bz
            else:
                self.bzs = -self.bz
                self.bzn = 0.0
            self.fsw = self.bperp*np.sqrt(np.sin(0.5*theta*self.__dtor))

        return

    def get_pressure(self,x,y):
        rho = np.sqrt(x*x+y*y)/self.__rnorm
        rm1 = rho-1.0
        phi = np.arctan2(y,x)

        P = (self.__aP[1]*rho**self.__aP[6]+
             self.__aP[2]*rho**self.__aP[7]*self.psw**self.__aP[11]+
             self.__aP[3]*rho**self.__aP[8]*self.fsw**self.__aP[12]+
             np.sin(phi)**2*
             (
                self.__aP[4]*self.psw**self.__aP[13]*np.exp(-self.__aP[9]*rho)+
                self.__aP[5]*self.fsw**self.__aP[14]*np.exp(-self.__aP[1]*rho)
             )
            )

        return P

    def get_density(self,x,y):
        rho = np.sqrt(x*x+y*y)/self.__rnorm
        rm1 = rho-1.0
        phi = np.arctan2(y,x)

        N = (rho**self.__aN[8]*
                (
                    self.__aN[1]+
                    self.__aN[2]*self.nsw**self.__aN[10]+
                    self.__aN[3]*self.bzn+
                    self.__aN[4]*self.vsw*self.bzs
                )+
                rho**self.__aN[9]*np.sin(phi)**2*
                (
                    self.__aN[5]*self.nsw**self.__aN[11]+
                    self.__aN[6]*self.bzn+
                    self.__aN[7]*self.vsw*self.bzs
                )
            )

        return N


    def get_temperature(self,x,y):
        rho = np.sqrt(x*x+y*y)/self.__rnorm
        rm1 = rho-1.0
        phi = np.arctan2(y,x)

        T = (self.__aT[1]*self.vsw+self.__aT[2]*self.bzn+self.__aT[3]*self.bzs+
             self.__aT[4]*np.exp(-rm1*
             (
                 self.__aT[9]*self.vsw**self.__aT[15]+
                 self.__aT[10]*self.bzn+
                 self.__aT[11]*self.bzs
             ))+
             np.sin(phi)**2*
             (
                self.__aT[5]*self.vsw+
                self.__aT[6]*self.bzn+
                self.__aT[7]*self.bzs+
                self.__aT[8]*np.exp(-rm1*
                (
                    self.__aT[12]*self.vsw**self.__aT[16]+
                    self.__aT[13]*self.bzn+
                    self.__aT[14]*self.bzs
                )
             ))
            )

        return T

    def calculate_moments(self,x,y):

        T = self.get_temperature(x,y)
        N = self.get_density(x,y)
        P = self.get_pressure(x,y)

        self.T = T
        self.N = N
        self.P = P

        self.moments = {'N':self.N,'P':self.P,'T':self.T}

        return self.moments

if __name__ == "__main__":

    import matplotlib.pyplot as pl

    xwant = np.linspace(-8,-48,41)
    ywant = np.linspace(-5,5,51)

    pmat = np.zeros([xwant.size,ywant.size])
    nmat = np.zeros_like(pmat)
    tmat = np.zeros_like(pmat)

    tm03 = plasma_sheet()
    tm03.set_parameters(vx=600,p=8.0)

    for i,x in enumerate(xwant):
        for j,y in enumerate(ywant):
            res=tm03.calculate_moments(x,y)
            pmat[i,j] = res['P']
            nmat[i,j] = res['N']
            tmat[i,j] = res['T']

    ax1 = pl.subplot2grid((3,1),(0,0))
    ax2 = pl.subplot2grid((3,1),(1,0))
    ax3 = pl.subplot2grid((3,1),(2,0))

    cp1=ax1.contourf(-xwant,ywant,(pmat.T),30,cmap=pl.get_cmap('nipy_spectral'))
    pl.colorbar(cp1,ax=ax1)
    cp2=ax2.contourf(-xwant,ywant,(nmat.T),30,cmap=pl.get_cmap('nipy_spectral'))
    pl.colorbar(cp2,ax=ax2)
    cp3=ax3.contourf(-xwant,ywant,(tmat.T),30,cmap=pl.get_cmap('nipy_spectral'))
    pl.colorbar(cp3,ax=ax3)

    pl.tight_layout()

    pl.show()