#!/usr/bin/env python
"""
SHIELDS-PTM run setup script for SEP access
"""
import argparse
import os
import subprocess
from contextlib import contextmanager
import numpy as np

from scripts import ptm_input


@contextmanager
def cd(newdir):
    """Context-managed 'change directory' """
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def setupGPS(opt, runid, nruns, verbose=False):
    # Get size of grid
    indir = os.path.join(opt.input_dir, 'ptm_data')
    # TODO: test for presence of directory and files, if not there, error
    xs = np.fromfile(os.path.join(indir, 'xgrid.bin'))
    ys = np.fromfile(os.path.join(indir, 'ygrid.bin'))
    zs = np.fromfile(os.path.join(indir, 'zgrid.bin'))

    # break energy space over number of runs...
    emin = 5.0*1e3  # 5MeV lower limit
    emax = 100.0*1e3  # 100MeV upper limit
    nenergy = 100
    energy_arr = np.linspace(emin, emax, nenergy)
    per_run = nenergy//nruns

    for runno in range(nruns):
        imin = per_run*runno
        imax = per_run*(runno+1)  # exclusive for use in ranges
        earr = {'min': energy_arr[imin],
                'max': energy_arr[imax-1],
                'n_e': len(energy_arr[imin:imax])
                }
        # TODO: presumably will need to do multiple run definitions here...
        #       then map them all to a GNU parallel launcher
        setup = ptm_input.ptm_input_creator(runid=runid+runno,
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
        setup.set_parameters(emin=earr['min'],  # Emin in keV
                             emax=earr['max'],  # Emax in keV
                             phi=-1.0,  # negative phi is randomly-seeded
                             nenergy=earr['n_e'],  # number of energies
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

        # Write required input files to run directory
        setup.create_input_files(os.path.join(opt.output_dir, 'ptm_input'))

        # Verbose output
        if verbose:
            setup.print_settings()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--runid', dest='runid', type=int, default=2,
                        help='Set starting run ID number. Default is 2.')
    parser.add_argument('-i', '--input', dest='input_dir',
                        default=None,
                        type=os.path.expanduser,
                        help='Path to input data files. Default is None. '
                        + 'If not provided, the script will assumes that the '
                        + 'preprocessing is done and in [output_dir]/ptm_data')
    parser.add_argument('-o', '--output', dest='output_dir',
                        default=os.path.abspath('.'),
                        type=os.path.expanduser,
                        help='Path to PTM run directory files. Default is "." '
                        + 'Params will go in [output_dir]/ptm_input. '
                        + 'Pre-processed data files go in [output_dir]/ptm_data.')
    parser.add_argument('--tec', dest='tec', action='store_true',
                        help='Set flag if SWMF input needs to be processed from '
                        + 'Tecplot format. Default is IDL format.')
    parser.add_argument('-t', '--timestep', dest='timestep', type=float, default=300.0)
    parser.add_argument('-p', '--position', dest='start_pos', nargs=3,
                        type=float)

    opt = parser.parse_args()

    # make sure run directory is built
    for dirm in ['ptm_data', 'ptm_input', 'ptm_output']:
        os.makedirs(os.path.join(opt.output_dir, dirm),
                    mode=0o775, exist_ok=True)

    # make symlink to ptm executable
    exe = os.path.abspath('ptm')
    link_to = os.path.join(opt.output_dir, 'ptm')
    if os.path.isfile(exe) and not os.path.isfile(link_to):
        os.symlink(exe, link_to)

    # if required, do preprocessing
    if opt.input_dir is None:
        opt.input_dir = opt.output_dir
    else:
        swmf_in = os.path.expandvars(os.path.expanduser(opt.input_dir))
        arglist = ["python", "ptm_preprocessing.py",
                   "-i", swmf_in,
                   "-o", os.path.join(opt.output_dir, 'ptm_data'),
                   "-t", '{0}'.format(opt.timestep),
                   ]
        if opt.tec:
            arglist.append("--tec")
        with cd('scripts'):
            print(os.path.abspath(os.curdir))
            subprocess.run(arglist)
            print('Done preprocessing')
        opt.input_dir = opt.output_dir

    # Do setup
    # TODO: break up into multiple runs to pass to GNU parallel
    nRuns = 3
    setupGPS(opt, opt.runid, nRuns)

    # TODO: have local run option as well as HPC run option (job script)
    if False:
        # run on non-cluster machine (e.g. desktop, scheme server)
        with cd(opt.output_dir): 
            subprocess.run(["ptm", "{0}".format(rundir)])
