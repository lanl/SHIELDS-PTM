"""
This script file can be used to interpolate data from an SMWF run. It is assumed
that the data is in gzipped tecplot-formatted ascii files.
"""

import glob
import os
import numpy as np
import ptm_read
import ptm_interpolate

# If you don't sort first, file order will be unpredictable
files = np.sort(glob.glob('*.gz'))

dtout = 300.0
times = dtout*np.arange(np.size(files))
times.tofile('ptm_data/tgrid.bin')

xwant = np.linspace(-15,10,101)
ywant = np.linspace(-10,10,81)
zwant = np.linspace(-10,10,81)

xwant.tofile('ptm_data/xgrid.bin')
ywant.tofile('ptm_data/ygrid.bin')
zwant.tofile('ptm_data/zgrid.bin')

for i, fname in enumerate(files):

    # The stat variable can be used to handle exceptions
    stat=os.system('gunzip '+fname)

    # Parse the newly unzipped file
    dat=ptm_read.read_swmf_tec_file(fname[:-3])

    # Rezip the file
    stat=os.system('gzip '+fname[:-3])

    # Interpolate the data
    res=ptm_interpolate.gauss_interp_EB(xwant,ywant,zwant,dat)

    res['bx'].tofile('ptm_data/bx3d_%4.4i.bin' % i)
    res['by'].tofile('ptm_data/by3d_%4.4i.bin' % i)
    res['bz'].tofile('ptm_data/bz3d_%4.4i.bin' % i)
    res['ex'].tofile('ptm_data/ex3d_%4.4i.bin' % i)
    res['ey'].tofile('ptm_data/ey3d_%4.4i.bin' % i)
    res['ez'].tofile('ptm_data/ez3d_%4.4i.bin' % i)