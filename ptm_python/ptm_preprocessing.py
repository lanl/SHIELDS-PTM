"""
This script file can be used to interpolate data from an SMWF run. It is assumed
that the data is in SWMF 3d idl format (binary OR ascii)
Binary is preferred due to the smaller size on disk
The ASCII tecplot format is also supported.
"""
import argparse
import glob
import os
import numpy as np
# local directory imports
import ptm_read
import ptm_interpolate

# Python 2 compatibility for checking string types
if not hasattr(__builtins__, "basestring"):
    basestring = (str, bytes)

def convertSWMF(files, dest_dir, dt=300, tec=False, extents='sep', static=False):
    """Converts input SWMF file to PTM-required binary format

    Parameters
    ==========
    files : list-like
        List of input filenames to process
    tec : boolean
        Keyword argument, select True for TecPlot format,
        False for IDL ASCII or binary
    extents : dict-like or string
        If dict-like input, extents can be manually given
        Requires 'xmin', 'xmax', 'xres', and similar for 'y' and
        'z'. Each value should be numeric. E.g.,
        extents = {'xmin': -15, 'xmax': 15, 'xres': 121,
                   'ymin': -15, 'ymax': 15, 'yres': 121,
                   'zmin': -15, 'zmax': 15, 'zres': 121}
        If string input, options are:
            'sep' : [-15, 15] in all dimensions, 121 divisions
            'geo' : [-10, 10] in (y,z), 81 divisions;
                    [-15, 10] in x, 101 divisions
    static : boolean
        Keyword argument, select True to duplicate a single timestep
    """
    times = dt*np.arange(np.size(files))
    times.tofile(os.path.join(dest_dir, 'tgrid.bin'))

    # Set resolution of data cube for PTM to trace through
    reqkeys = ['xmin', 'ymin', 'zmin',
               'xmax', 'ymax', 'zmax',
               'xres', 'yres', 'zres']
    missing = []
    if isinstance(extents, dict):
        for key in reqkeys:
            if key not in extents:
                missing.append(key)
        if missing:
            raise KeyError('convertSWMF: extents dict requires min, max and res '
                           + 'to be set for each of x, y, z\n'
                           + 'Missing: {0}'.format(missing))
        extdict = extents
    elif isinstance(extents, basestring):
        if extents.lower() == 'sep':
            extdict = {'xmin': -15, 'xmax': 15, 'xres': 121,
                       'ymin': -15, 'ymax': 15, 'yres': 121,
                       'zmin': -15, 'zmax': 15, 'zres': 121}
        elif extents.lower() == 'geo':
            extdict = {'xmin': -15, 'xmax': 10, 'xres': 101,
                       'ymin': -10, 'ymax': 10, 'yres': 81,
                       'zmin': -10, 'zmax': 10, 'zres': 81}
        else:
            raise ValueError('convertSWMF: extents must be provided as dict or name')
    else:
        raise ValueError('convertSWMF: extents must be provided as dict or name')
    xwant = np.linspace(extdict['xmin'], extdict['xmax'], extdict['xres'])
    ywant = np.linspace(extdict['ymin'], extdict['ymax'], extdict['yres'])
    zwant = np.linspace(extdict['zmin'], extdict['zmax'], extdict['zres'])

    # tofile is a numpy method that dumps binary (which the fortran will read)
    xwant.tofile(os.path.join(dest_dir, 'xgrid.bin'))
    ywant.tofile(os.path.join(dest_dir, 'ygrid.bin'))
    zwant.tofile(os.path.join(dest_dir, 'zgrid.bin'))

    if static:
        files = [files[0], files[0]]

    for i, fname in enumerate(files):
        # The stat variable can be used to handle exceptions
        #stat=os.system('gunzip '+fname)

        # Parse the newly unzipped file
        #dat=ptm_read.read_swmf_tec_file(fname[:-3])
        doRead = (i == 0) or (i > 0 and not static)
        if doRead:
            print('start reading', fname)
            if tec:
                # TODO: Also reinstate gzip handling
                dat = ptm_read.read_swmf_tec_file(fname)
            else:
                dat = ptm_read.read_swmf_idl_file(fname)
            print('done reading file', fname)

            # Interpolate the data
            res = ptm_interpolate.gauss_interp_EB(xwant, ywant, zwant, dat)
            print('done interpolating file', fname)

        i1 = i+1
        res['Bx'].tofile(os.path.join(dest_dir, 'bx3d_{0:04d}.bin'.format(i1)))
        res['By'].tofile(os.path.join(dest_dir, 'by3d_{0:04d}.bin'.format(i1)))
        res['Bz'].tofile(os.path.join(dest_dir, 'bz3d_{0:04d}.bin'.format(i1)))
        res['Ex'].tofile(os.path.join(dest_dir, 'ex3d_{0:04d}.bin'.format(i1)))
        res['Ey'].tofile(os.path.join(dest_dir, 'ey3d_{0:04d}.bin'.format(i1)))
        res['Ez'].tofile(os.path.join(dest_dir, 'ez3d_{0:04d}.bin'.format(i1)))

if __name__ == '__main__':
    # Define command-line option parser
    curr = os.path.abspath(os.path.curdir)
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_dir', dest='input_dir',
                        default=os.path.join(curr, '..', 'ptm_data'))
    parser.add_argument('-o', '--output_dir', dest='output_dir',
                        default=os.path.join(curr, '..', 'ptm_data'))
    parser.add_argument('-t', '--timestep', dest='timestep',
                        default=300, type=float,
                        help='Time between E/B snapshots in seconds. Default = 300')
    parser.add_argument('--tec', dest='tec', action='store_true',
                        help='Flag to process as TecPlot format. Default is IDL format.')
    parser.add_argument('--geo', dest='gridtype', action='store_true',
                        help='Flag to process onto original "geo" grid. Default is new "sep" grid.')
    parser.add_argument('fname', nargs='*', help='Filenames to convert to PTM format. '
                        + 'Any number can be provided. Wildcards assumed to be expanded by shell.\n'
                        + 'Default is to use *.out from input_dir.')
    opt = parser.parse_args()
    grid = 'geo' if opt.gridtype else 'sep'
    if not opt.fname:
        opt.fname = glob.glob(os.path.join(opt.input_dir, '*.out'))
    if not os.path.isdir(opt.output_dir):
        os.makedirs(opt.output_dir)
    # If you don't sort first, file order will be unpredictable
    files = np.sort(opt.fname)
    # If only one timestep, make sure we have the file input twice...
    # Allows in-time interpolation in a static field
    static = True if len(files) == 1 else False
    convertSWMF(files, opt.output_dir, tec=opt.tec, extents=grid, static=static)
