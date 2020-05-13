#!/usr/bin/env python
"""
SHIELDS-PTM run setup script for SEP access
"""
import argparse
import os
import numpy as np

from scripts import ptm_input

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input_dir',
                        default=os.path.abspath('ptm_data'),
                        help='Path to input data files. Default is ptm_data.')
    parser.add_argument('-p', '--position', dest='start_pos', nargs=3,
                        type=float)

    opt = parser.parse_args()
    print(opt.start_pos)
    # Get size of grid
    # TODO: test for presence of directory and files, if not there, error
    xs = np.fromfile(os.path.join(opt.input_dir, 'xgrid.bin'))
    ys = np.fromfile(os.path.join(opt.input_dir, 'ygrid.bin'))
    zs = np.fromfile(os.path.join(opt.input_dir, 'zgrid.bin'))

    # TODO: presumably will need to do multiple run definitions here...
    #       then map them all to a GNU parallel launcher
    setup = ptm_input.ptm_input_creator(runid=2,
                                        idensity=1,  # Point source
                                        )

    # Ensure domain resolution is correctly set
    setup.set_parameters(ndim=3,
                         nx=xs.shape[0],
                         ny=ys.shape[0],
                         nz=zs.shape[0])

    # Change to flux map mode
    setup.set_parameters(idist=3)  # TODO: can I get particle 
                                   # trajectories in the usual
                                   # format with this mode?
    # How do I combine multiple "runs"? For efficiency I should
    # parcel the runs up, probably using smaller ranges in
    # energy...
    setup.set_parameters(emin=10000.0,  # Emin in keV [10MeV]
                         emax=100000.0,  # Emax in keV [100MeV]
                         phi=-1.0,  # negative phi is randomly-seeded
                         nenergy=10,  # number of energies
                         npitch=18,  # number of pitch angles
                         xsource=-20.0
                         )

    # Set proton tracing params, incl. full orbit
    setup.set_parameters(iswitch=1,  # full orbit
                         itrace=-1,  # backwards in time
                         charge=1.0, mass=1837.0,  # protons
                         dtout=0.1, # 0.1s output written
                         thi=300.0, # max time of 300s
                         itraj=1,  # 1=write trajectories when in fluxmap mode
                         )

    # Define point source for tracing
    setup.set_parameters(x0=opt.start_pos[0],
                         y0=opt.start_pos[1],
                         z0=opt.start_pos[2])

    setup.create_input_files('ptm_input')
    setup.print_settings()
