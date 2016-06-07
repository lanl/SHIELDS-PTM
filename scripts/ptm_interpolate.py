#!/usr/bin/python

"""
+
PTM_INTERPOLATE

Interpolation of SWMF data files to create input data for the ptm particle tracing routines.
+

This file is intended for use as a command line script although there are individual routines within
that can be used interactively.

COMMAND LINE USAGE

./ptm_interpolate 3d_mhd__00100000_0024350.out 42

This will read the fields from the file "3d_mhd__..." and will write out six files:

bx3d_0042.bin
by3d_0042.bin
bz3d_0042.bin
ex3d_0042.bin
ey3d_0042.bin
ez3d_0042.bin

Additionally, the file 'tgrid.bin' will be either updated to include the latest timstep or created if 
it does not already exist (in the example above case, 36000.0 would be added to the file).

CONTAINS

  gauss_interp_EB   Interpolate fields using a Gaussian Inverse Distance Weighting (IDW) method
  rbf_interp_EB     Interpolate fields using a Radial Basis Function (RBF) method

See individiual functions for descriptions and usage.

Jesse Woodroffe
last updated 6/2/2016

"""

from numpy import meshgrid,tile,mean,exp,sum,reshape,einsum,fromfile,c_,newaxis,size,zeros_like,shape
from spacepy.pybats import IdlFile
from scipy.spatial import cKDTree
from scipy.interpolate import Rbf
from sys import argv
from os.path import isfile

def gauss_interp_EB(xwant,ywant,zwant,fdict,smoothingDegree=0.5,numNeighbors=24):
  """
  Three-dimensional Gaussian weighting interpolator. This is a form of inverse distance weighting (IDW)
  that looks a lot like a 3D kernal density estimator. This particular implementation has been designed
  to interpolate data (B and E=-u x B) from the scattered SWMF mesh onto a uniform grid.
   
  Inputs:
    xwant   Array of x-positions where values are desired
    ywant   Array of y-positions where values are desired
    zwant   Array of z-positions where values are desired
    fdict   Dictionary containing data from an SWMF object created by IdlFile
    
  Optional:
    smoothingDegree   Scalar value used to determine the strength of contributions from distant points
    numNeighbors      Integer value setting the number of neighbors to use in calculation
    
  Jesse Woodroffe 5/20/2016
  """
  
  mhdFields=['bx','by','bz','ux','uy','uz']
  
  # Check that the data file is three-dimensional
  try:
    myTree = cKDTree(zip(fdict['x'],fdict['y'],fdict['z']))
  except KeyError:
    missingKeys=''
    if(not fdict.has_key('x')): missingKeys+='x '
    if(not fdict.has_key('y')): missingKeys+='y '
    if(not fdict.has_key('z')): 
      missingKeys+='z'
      if(missingKeys=='z'): print "\nBased on missing keys, SWMF file may contain only 2D data\n"
    raise Exception('Error in rbf_interp_EB: fdict does not have expected keys: {'+missingKeys.strip()+'}')
    
  # Check that data file contains the correct components  

  for key in mhdFields:
    if(not fdict.has_key(key)): raise Exception('Error in rbf_interp_EB: requested key {'+key+'} was not present in fdict')
    
  Y,X,Z=meshgrid(ywant,xwant,zwant)
  dists,dexes=myTree.query(c_[X.ravel(),Y.ravel(),Z.ravel()],numNeighbors)
  epss=smoothingDegree*mean(dists,axis=1)
  scal=dists*tile(1.0/epss[:,newaxis],size(dists,1))
  sbfs=exp(-scal*scal)
  sbms=1.0/sum(sbfs,axis=1)
    
  temp={}
    
  for key in mhdFields:
    temp[key]=reshape(einsum('...j,...j',sbfs,fdict[key][dexes])*sbms,X.shape)

  res={}
  res['bx']=temp['bx'].copy()
  res['by']=temp['by'].copy()
  res['bz']=temp['bz'].copy()
  res['ex']=-(temp['uy']*temp['bz']-temp['uz']*temp['by'])
  res['ey']=-(temp['uz']*temp['bx']-temp['ux']*temp['bz'])
  res['ez']=-(temp['ux']*temp['by']-temp['uy']*temp['bx'])

  return res

def rbf_interp_EB(xwant,ywant,zwant,fdict,numNeighbors=48,smoothingDegree=0.05,basis='multiquadric'):
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
    fdict   Dictionary containing data from an SWMF object created by IdlFile from spacepy.pybats

  Optional:
    numNeighbors      Integer value setting the number of neighbors to use in calculation, default is 48
    basis             String specifying the Rbf interpolation kernel to be used, default is multiquadric
    smoothingDegree   Float that determines how closesly input values need to be replicated (0.0=exactly)
  
  Jesse Woodroffe
  6/2/2016
  """
  
  goodBases=['multiquadric','inverse','thin_plate','gaussian','linear','cubic','quintic']
  mhdFields=['bx','by','bz','ex','ey','ez']
  
  # Check that user is passing three-dimensional data
  try:
    myTree = cKDTree(zip(fdict['x'],fdict['y'],fdict['z']))
  except KeyError:
    missingKeys=''
    if(not fdict.has_key('x')): missingKeys+='x '
    if(not fdict.has_key('y')): missingKeys+='y '
    if(not fdict.has_key('z')): 
      missingKeys+='z'
      if(missingKeys=='z'): print "\nBased on missing keys, SWMF file may contain only 2D data\n"
    raise Exception('Error in rbf_interp_EB: fdict does not have expected keys: {'+missingKeys.strip()+'}')

  # Check that requested basis is supported  
  if(not basis in goodBases):
    raise Exception('Error in rbf_interp_EB: requested basis ('+basis+') not supported')

  # Check that data file contains the correct components  
  for key in mhdFields:
    if(not fdict.has_key(key)): raise Exception('Error in rbf_interp_EB: requested key {'+key+'} was not present in fdict')
  
  myTree = cKDTree(zip(fdict['x'],fdict['y'],fdict['z']))
  Y,X,Z=meshgrid(ywant,xwant,zwant)
  XR,YR,ZR=X.ravel(),Y.ravel(),Z.ravel()
  dists,dexes=myTree.query(c_[XR,YR,ZR],numNeighbors)
  
  temp={}
  for key in goodFields:
    temp[key]=zeros_like(ZR)  
    for i in xrange(size(ZR)):
      myRBF = Rbf(fdict['x'][dexes[i,:]],fdict['y'][dexes[i,:]],fdict['z'][dexes[i,:]],fdict[key][dexes[i,:]],smooth=smoothingDegree,function=basis)
      temp[key][i] = myRBF(XR[i],YR[i],ZR[i])

  res={}
  res['bx']=temp['bx'].reshape(Z.shape)
  res['by']=temp['by'].reshape(Z.shape)
  res['bz']=temp['bz'].reshape(Z.shape)
  res['ex']=-(temp['uy']*temp['bz']-temp['uz']*temp['by']).reshape(Z.shape)
  res['ey']=-(temp['uz']*temp['bx']-temp['ux']*temp['bz']).reshape(Z.shape)
  res['ez']=-(temp['ux']*temp['by']-temp['uy']*temp['bx']).reshape(Z.shape)

  return res

if __name__ == "__main__":

  """
  If called as a script with arguments <filename> <id number>, interpolates E and B and creates
  output files. This script should be called once for every new SWMF data file that is created.
  
  Jesse Woodroffe 5/20/2016
  """

  # Determine time of data output by parsing the input file name
  filename=argv[1]
  tstring=filename.split('_')[-2][1:]
  dd=86400*float(tstring[:2])
  hh=3600*float(tstring[2:4])
  mm=60.0*float(tstring[4:6])
  ss=float(tstring[6:])
  time=dd+hh+mm+ss

  try:
    swmfdata = IdlFile(filename)
    xgrid = fromfile('xgrid.bin')
    ygrid = fromfile('ygrid.bin')
    zgrid = fromfile('zgrid.bin')
  except FileNotFoundError:
    if(not isfile('xgrid.bin')): print 'xgrid.bin not found'
    if(not isfile('ygrid.bin')): print 'ygrid.bin not found'
    if(not isfile('zgrid.bin')): print 'zgrid.bin not found'
    if(not isfile(filename)): print filename+' not found'
    raise
  
  # We could use either Gauss or Rbf interpolator, Gauss is faster so we'll go with it.   
  fields=gauss_interp_EB(xgrid,ygrid,zgrid,swmfdata)
 
  for key in ['bx','by','bz','ex','ey','ez']:
    fields[key].tofile(key+'3d_%4.4i.bin' % argv[2])
    
  if(isfile('tgrid.bin')):
    # Read in the tgrid array and add the new time value to
    tgrid = fromfile('tgrid.bin')
    tgrid.append(time)
  else:
    # Create the tgrid array
    tgrid=r_[time]
    
  tgrid.tofile('tgrid.bin')
