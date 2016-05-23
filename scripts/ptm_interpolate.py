#!/usr/bin/python

"""
Example use

./ptm_interpolate 3d_mhd__00100000_0024350.out 42

This will read the fields from the file "3d_mhd__..." and will write out six files:

bx3d_0042.bin
by3d_0042.bin
bz3d_0042.bin
ex3d_0042.bin
ey3d_0042.bin
ez3d_0042.bin

In addition, the file 'tgrid.bin' will be updated to include the latest timstep.

"""

from numpy import meshgrid,tile,mean,exp,sum,reshape,einsum,fromfile
from spacepy.pybats import IdlFile
from scipy.spatial import cKDTree
from sys import argv
from os.path import isfile

def gauss_interp_3xN(xwant,ywant,zwant,fdict,smoothingDegree=0.5,numNeighbors=24):
    """
    Three-dimensional Gaussian weighting interpolator for N-component irregularly gridded data. This is a form
    of inverse distance weighting that looks a lot like a 3D kernal density estimator.
    
    Jesse Woodroffe 5/20/2016
    """
    
    myTree = cKDTree(zip(fdict[0]['x'],fdict[0]['y'],fdict[0]['z']))
    
    Y,X,Z=meshgrid(ywant,xwant,zwant)
    dists,dexes=dataTree.query(c_[X.ravel(),Y.ravel(),Z.ravel()],numNeighbors)
    epss=smoothingDegree*mean(dists,axis=1)
    scal=dists*tile(1.0/epss[:,newaxis],size(dists,1))
    sbfs=exp(-scal*scal)
    sbms=1.0/sum(sbfs,axis=1)
    
    temp={}
    
    for key in ['bx','by','bz','ux','uy','uz']:
      temp[key]=reshape(einsum('...j,...j',sbfs,fdict[key][dexes])*sbms,X.shape)

    res={}
    res['bx']=temp['bx'].copy()
    res['by']=temp['by'].copy()
    res['bz']=temp['bz'].copy()
    res['ex']=-(temp['uy']*temp['bz']-temp['uz']*temp['by'])
    res['ey']=-(temp['uz']*temp['bx']-temp['ux']*temp['bz'])
    res['ez']=-(temp['ux']*temp['by']-temp['uy']*temp['bx'])

    return res
    
if __name__ == "__main__":

  """
  If called as a script with arguments <filename> <id number>, interpolates E and B and creates
  output files. This script should be called once for every new SWMF data file that is created.
  
  Jesse Woodroffe 5/20/2016
  """

  swmfdata = IdlFile(argv[1])
  xgrid = fromfile('xgrid.bin')
  ygrid = fromfile('ygrid.bin')
  zgrid = fromfile('zgrid.bin')
  
  # We could use a different interpolator if one becomes available    
  fields=gauss_interp_3xN(xgrid,ygrid,zgrid,swmfdata)
 
  for key in ['bx','by','bz','ex','ey','ez']:
    fields[key].tofile(key+'3d_%4.4i.bin' % argv[2])
    
  # Update the time file to include this most recent step
  if(isfile('tgrid.bin')):
    tgrid = fromfile('tgrid.bin')
    tgrid.append(swmfdata['time'])
  else:
    tgrid=r_[swmfdata['time']]
    
  tgrid.tofile('tgrid.bin')