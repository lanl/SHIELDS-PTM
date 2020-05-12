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
import ptm_read
import ptm_interpolate

def convertSWMF(files, dest_dir, tec=False):
    dtout = 3600.
    times = dtout*np.arange(np.size(files))
    times.tofile(os.path.join(dest_dir, 'tgrid.bin'))

    # Set resolution of data cube for PTM to trace through
    xwant = np.linspace(-15,15,75)
    ywant = np.linspace(-15,15,75)
    zwant = np.linspace(-15,15,75)

    # tofile is a numpy method that dumps binary (which the fortran will read)
    xwant.tofile(os.path.join(dest_dir, 'xgrid.bin'))
    ywant.tofile(os.path.join(dest_dir, 'ygrid.bin'))
    zwant.tofile(os.path.join(dest_dir, 'zgrid.bin'))

    for i, fname in enumerate(files):
        # The stat variable can be used to handle exceptions
        #stat=os.system('gunzip '+fname)

        # Parse the newly unzipped file
        #dat=ptm_read.read_swmf_tec_file(fname[:-3])
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
    parser.add_argument('--tec', dest='tec', action='store_true',
                        help='Flag to process as TecPlot format. Default is IDL format.')
    parser.add_argument('fname', nargs='*', help='Filenames to convert to PTM format. '
                        + 'Any number can be provided. Wildcards assumed to be expanded by shell.\n'
                        + 'Default is to use *.out from input_dir.')
    opt = parser.parse_args()
    if not opt.fname:
        opt.fname = glob.glob(os.path.join(opt.input_dir, '*.out'))
    if not os.path.isdir(opt.output_dir):
        os.makedirs(opt.output_dir)
    # If you don't sort first, file order will be unpredictable
    files = np.sort(opt.fname)
    if len(files) == 1:
        # If only one timestep, make sure we have the file input twice...
        # Allows in-time interpolation in a static field
        files = [files[0], files[0]]
    convertSWMF(files, opt.output_dir, tec=opt.tec)
