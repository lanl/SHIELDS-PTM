#!/usr/bin/env python

"""
PTM_INTERPOLATE

Interpolation routines for 3D scattered mesh data that can be used for preparing SWMF output
to be used by the PTM code.

CONTAINS

  gauss_interp_EB   Interpolate fields using a Gaussian Inverse Distance Weighting (IDW) method
  rbf_interp_EB     Interpolate fields using a Radial Basis Function (RBF) method

See individual functions for descriptions and usage.

Jesse Woodroffe, Steven Morley
last updated 5/11/2020
"""

import numpy as np
from scipy.spatial import cKDTree
from scipy.interpolate import Rbf

__mhdFields = ['Bx','By','Bz','Ux','Uy','Uz']


def makeTree(fdict, calling='unknown func'):
    # Check that the data file is three-dimensional
    try:
        inarray = np.empty((len(fdict['x']), 3))
        inarray[:, 0] = fdict['x'].squeeze()
        inarray[:, 1] = fdict['y'].squeeze()
        inarray[:, 2] = fdict['z'].squeeze()
        myTree = cKDTree(inarray)
    except KeyError:
        missingKeys=''
        if (not fdict.has_key('x')):
            missingKeys += 'x '
        if (not fdict.has_key('y')):
            missingKeys+='y '
        if (not fdict.has_key('z')):
            missingKeys+='z'
        raise Exception('Error in {0}: '.format(calling)
                        + 'fdict does not have expected keys: {'
                        + missingKeys.strip()+ '}')

    # Check that data file contains the correct components
    for key in __mhdFields:
        if (key not in fdict):
            raise Exception('Error in {0}: requested key'.format(calling)
                            + ' {' + key + '} was not present in fdict')

    return myTree


def gauss_interp(xwant, ywant, zwant, fdict, key, smoothingDegree=0.5, numNeighbors=24):
    """
    Three-dimensional Gaussian weighting interpolator. This is a form of inverse distance weighting (IDW)
    that looks a lot like a 3D kernal density estimator.

    Inputs:
        xwant   Array of x-positions where values are desired
        ywant   Array of y-positions where values are desired
        zwant   Array of z-positions where values are desired
        fdict   Dictionary containing data from an SWMF run

    Optional:
        smoothingDegree   Scalar value used to determine the strength of contributions from distant points
        numNeighbors      Integer value setting the number of neighbors to use in calculation
    """
    myTree = makeTree(fdict, calling='gauss_interp')
    Y, X, Z = np.meshgrid(ywant, xwant, zwant)
    dists, dexes = myTree.query(np.c_[X.ravel(), Y.ravel(), Z.ravel()], numNeighbors)
    epss = smoothingDegree*np.mean(dists, axis=1)
    scal = dists*np.tile(1.0/epss[:, np.newaxis], np.size(dists, 1))
    sbfs = np.exp(-scal*scal)
    sbms = 1.0/np.sum(sbfs, axis=1)

    res = np.reshape(np.einsum('...j,...j', sbfs, fdict[key].squeeze()[dexes])*sbms, X.shape)

    return res


def gauss_interp_EB(xwant, ywant, zwant, fdict, smoothingDegree=0.5, numNeighbors=24):
    """
    Three-dimensional Gaussian weighting interpolator. This is a form of inverse distance weighting (IDW)
    that looks a lot like a 3D kernal density estimator. This particular implementation has been designed
    to interpolate data (B and E=-u x B) from the scattered SWMF mesh onto a uniform grid.

    Inputs:
        xwant   Array of x-positions where values are desired
        ywant   Array of y-positions where values are desired
        zwant   Array of z-positions where values are desired
        fdict   Dictionary containing data from an SWMF run

    Optional:
        smoothingDegree   Scalar value used to determine the strength of contributions from distant points
        numNeighbors      Integer value setting the number of neighbors to use in calculation
    """
    myTree = makeTree(fdict, calling='gauss_interp_EB')
    Y, X, Z = np.meshgrid(ywant, xwant, zwant)
    dists, dexes = myTree.query(np.c_[X.ravel(), Y.ravel(), Z.ravel()], numNeighbors)
    epss = smoothingDegree*np.mean(dists, axis=1)
    scal = dists*np.tile(1.0/epss[:, np.newaxis], np.size(dists, 1))
    sbfs = np.exp(-scal*scal)
    sbms = 1.0/np.sum(sbfs,axis=1)

    temp = {}

    for key in __mhdFields:
        temp[key] = np.reshape(np.einsum('...j,...j', sbfs, fdict[key].squeeze()[dexes])*sbms, X.shape)

    res = {}
    res['Bx'] = temp['Bx'].copy()
    res['By'] = temp['By'].copy()
    res['Bz'] = temp['Bz'].copy()
    res['Ex'] = -(temp['Uy']*temp['Bz']-temp['Uz']*temp['By'])
    res['Ey'] = -(temp['Uz']*temp['Bx']-temp['Ux']*temp['Bz'])
    res['Ez'] = -(temp['Ux']*temp['By']-temp['Uy']*temp['Bx'])

    return res


def rbf_interp_EB(xwant, ywant, zwant, fdict, numNeighbors=48, smoothingDegree=0.05, basis='multiquadric'):
    """
    Perform radial basis function interpolation of an irregularly-gridded data set using
    a K-nearest neighbors approach. This method is considerably  slower (approx 12x) than the IDW
    approach in gauss_interp_EB, owing to its need to solve a linear system for each point
    being interpolated. Consequently, this cannot be recommended for general use, but it is
    provided for comparison's sake. An additional issue is that Rbf methods are sensitive to
    resolution changes, producing a "ringing" effect near these boundaries. To this point, the
    basis that seems to deal best with these issues is the multiquadric, so this is default.

    Inputs:
        xwant   Array of x-positions where values are desired
        ywant   Array of y-positions where values are desired
        zwant   Array of z-positions where values are desired
        fdict   Dictionary containing data from an SWMF run

    Optional:
        numNeighbors      Integer value setting the number of neighbors to use in calculation, default is 48
        basis             String specifying the Rbf interpolation kernel to be used, default is multiquadric
        smoothingDegree   Float that determines how closesly input values need to be replicated (0.0=exactly)
    """

    goodBases=['multiquadric', 'inverse', 'thin_plate', 'gaussian',
               'linear', 'cubic', 'quintic']

    myTree = makeTree(fdict, calling='rbf_interp_EB')
    # Check that requested basis is supported
    if (not basis in goodBases):
        raise Exception('Error in rbf_interp_EB: requested basis ('
                        + basis + ') not supported')

    Y, X, Z = np.meshgrid(ywant, xwant, zwant)
    XR, YR, ZR = X.ravel(), Y.ravel(), Z.ravel()
    dists, dexes = myTree.query(c_[XR, YR, ZR], numNeighbors)

    temp={}
    for key in __mhdFields:
        temp[key] = np.zeros_like(ZR)
        for i in range(size(ZR)):
            myRBF = Rbf(fdict['x'][dexes[i, :]], fdict['y'][dexes[i, :]], fdict['z'][dexes[i, :]],
                        fdict[key][dexes[i, :]], smooth=smoothingDegree, function=basis)
            temp[key][i] = myRBF(XR[i], YR[i], ZR[i])

    res={}
    res['Bx'] = temp['Bx'].reshape(Z.shape)
    res['By'] = temp['By'].reshape(Z.shape)
    res['Bz'] = temp['Bz'].reshape(Z.shape)
    res['Ex'] = -(temp['Uy']*temp['Bz']-temp['Uz']*temp['By']).reshape(Z.shape)
    res['Ey'] = -(temp['Uz']*temp['Bx']-temp['Ux']*temp['Bz']).reshape(Z.shape)
    res['Ez'] = -(temp['Ux']*temp['By']-temp['Uy']*temp['Bx']).reshape(Z.shape)

    return res
