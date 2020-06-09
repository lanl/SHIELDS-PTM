#!/usr/bin/env python
"""
SHIELDS-PTM run setup script for SEP access
"""
import argparse
import os
import subprocess
from contextlib import contextmanager
import numpy as np

from ptm_python import ptm_input


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

    # break energy & alpha space over number of runs...
    emin = 5.0*1e3  # 5MeV lower limit
    emax = 800.0*1e3  # 800MeV upper limit
    nenergy = 300
    pamin = 0.1
    pamin_r = np.deg2rad(pamin)
    pamax = 179.9
    pamax_r = np.deg2rad(pamax)
    nalpha = 66
    energy_arr = np.logspace(np.log10(emin), np.log10(emax), num=nenergy)
    per_run = nenergy//nruns
    # pitch angles evenly spaced around unit sphere in velocity space
    alpha_arr = np.rad2deg(np.arccos(np.linspace(np.cos(pamin_r),
                                                 np.cos(pamax_r),
                                                 num=nalpha)))
    per_run_alpha = nalpha//nruns
    # now reorder energy and pitch angle arrays randomly
    np.random.shuffle(energy_arr)
    np.random.shuffle(alpha_arr)

    for runno in range(nruns):
        imin = per_run*runno
        imin_a = per_run_alpha*runno
        imax = per_run*(runno+1)  # exclusive for use in ranges
        imax_a = per_run_alpha*(runno+1)
        earr = {'n_e': per_run if runno != nruns-1 else len(energy_arr[imin:]),
                'n_a': per_run_alpha if runno != nruns-1 else len(alpha_arr[imin_a:]),
                }
        # build energy and pitch angle arrays and write to file
        runenergies = np.array(energy_arr[imin:imax])
        runalphas = np.array(alpha_arr[imin_a:imax_a])
        runenergies.tofile(os.path.join(indir, 'energies_{:04d}.bin'.format(runid+runno)))
        runalphas.tofile(os.path.join(indir, 'pitchangles_{:04d}.bin'.format(runid+runno)))

        # start input file setup
        setup = ptm_input.ptm_input_creator(runid=runid+runno,
                                            idensity=1,  # Point source
                                            )

        # Ensure domain resolution is correctly set
        setup.set_parameters(ndim=3,
                             nx=xs.shape[0],
                             ny=ys.shape[0],
                             nz=zs.shape[0])

        # Change to user-defined flux map mode
        setup.set_parameters(idist=4)
        # Break up into multiple "runs"? For efficiency I should
        # parcel the runs up, probably using smaller ranges in
        # energy...
        setup.set_parameters(phi=-1.0,  # negative phi is randomly-seeded
                             nenergy=earr['n_e'],  # number of energies
                             npitch=earr['n_a'],  # number of pitch angles
                             xsource=-20.0
                             )

        # Set proton tracing params, incl. full orbit
        setup.set_parameters(iswitch=1,  # full orbit
                             itrace=-1,  # backwards in time
                             charge=1.0, mass=1837.0,  # protons
                             dtout=0.1, # 0.1s output written
                             thi=opt.starttime, # max time ("start" for backwards trace)
                             tlo=opt.starttime-300,  # min time ("end" for backwards trace)
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


def setupElec(opt, runid, nruns, verbose=False):
    # Get size of grid
    indir = os.path.join(opt.input_dir, 'ptm_data')
    # TODO: test for presence of directory and files, if not there, error
    xs = np.fromfile(os.path.join(indir, 'xgrid.bin'))
    ys = np.fromfile(os.path.join(indir, 'ygrid.bin'))
    zs = np.fromfile(os.path.join(indir, 'zgrid.bin'))

    # break energy space over number of runs...
    emin = 50.0*1e-3  # 50eV lower limit
    emax = 200.0*1e-3  # 200eV upper limit
    nenergy = 30
    energy_arr = np.linspace(emin, emax, nenergy)
    per_run = nenergy//nruns

    for runno in range(nruns):
        imin = per_run*runno
        imax = per_run*(runno+1)  # exclusive for use in ranges
        earr = {'min': energy_arr[imin],
                'max': energy_arr[imax-1],
                'n_e': per_run if runno != nruns-1 else len(energy_arr[imin:])
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
                             npitch=9,  # number of pitch angles
                             pamin=0,
                             pamax=90,
                             xsource=-12.0
                             )

        # Set tracing params
        setup.set_parameters(iswitch=0,  # switch between GC and full orbit
                             itrace=-1,  # backwards in time
                             charge=-1.0, mass=1.0,  # electrons
                             dtout=1.0, # 0.1s output written
                             thi=600.0, # max time of 600s
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


def writeJobScript(options, nRuns, cluster=False, c_kwargs={}):
    runst = options.runid
    runen = options.runid + nRuns
    defaults = {'hh': 1,
                'mm': 0,
                'nodes': 3
                }
    if cluster:
        # Write a job script for SLURM
        for key in defaults:
            if key not in c_kwargs:
                c_kwargs[key] = defaults[key]
        headers = ['#!/usr/bin/env bash',
                   '#SBATCH --time={0:02d}:{1:02d}:00',
                   '#SBATCH --nodes={0}',
                   '#SBATCH --no-requeue',
                   '#SBATCH --job-name=ptm_parallel',
                   '#SBATCH -o slurm%j.out',
                   '#SBATCH -e slurm%j.err',
                   '#SBATCH --qos=standard',
                   '#SBATCH --account=w19_gic',
                   '#SBATCH --mail-user={0}@lanl.gov',
                   '#SBATCH --mail-type=FAIL',
                   ]
        body = ['module purge',
                'module load gcc',
                '',
                '# Export list of host (node) names to a file',
                'scontrol show hostnames > ./nodelist_${SLURM_JOB_ID}',
                '# record env to pass to remotes. Important for, e.g., OMP_NUM_THREADS',
                'export OMP_NUM_THREADS=16',
                '# parallel --record-env',
                ]
        with open('job_script.sh', 'w') as fh:
            for line in headers:
                if '--time' in line:
                    line.format(c_kwargs['hh'], c_kwargs['mm'])
                if '--nodes' in line:
                    line.format(c_kwargs['nodes'])
                if '--mail-user' in line:
                    line.format(os.environ['USER'])
                fh.write(line + '\n')
            fh.write('\n')
            for line in body:
                fh.write(line + '\n')
            pcmd = ' '.join(['parallel --jobs 2 --sshloginfile ./nodelist_${SLURM_JOB_ID}',
                             '--workdir=$PWD --env OMP_NUM_THREADS ./ptm :::',
                             '{' + '{0}..{1}'.format(runst, runen) + '}'
                             ])
            fh.write('\n' + pcmd + '\n')
        print('Submit jobs using `sbatch job_script.sh` from the run directory')
    else:
        with open('job_script.sh', 'w') as fh:
            fh.write('#!/usr/bin/env bash\n')
            cmds = ['./ptm {0}'.format(nn) for nn in range(runst, runen)]
            cmdstr = ' && '.join(cmds)
            fh.write(cmdstr + '\n')
        print('Run jobs using `source job_script.sh` from the run directory')


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
    parser.add_argument('-e', '--elec', dest='elec', action='store_true')
    parser.add_argument('-c', '--cluster', dest='cluster', action='store_true')
    parser.add_argument('-s', '--starttime', dest='starttime', type=float, default=300,
                        help='Start time for trajectory integration in integer seconds ' +
                        'since the start time of the fields input.')

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
        with cd('ptm_python'):
            print(os.path.abspath(os.curdir))
            subprocess.run(arglist)
            print('Done preprocessing')
        opt.input_dir = opt.output_dir

    # Do setup
    # TODO: break up into multiple runs to pass to GNU parallel
    if opt.elec:
        nRuns = 1
        setupElec(opt, opt.runid, nRuns)
    else:
        nRuns = 8
        setupGPS(opt, opt.runid, nRuns)

    # write job sscript: local run option as well as HPC run option
    if not opt.cluster:
        # run on non-cluster machine (e.g. desktop, scheme server)
        with cd(opt.output_dir):
            # local run should launch `./ptm n` for {opt.runid:nRuns}
            # TODO: write out to a script that will launch all jobs
            # subprocess.run(["ptm", "{0}".format(rundir)])
            writeJobScript(opt, nRuns, cluster=False)
    else:
        # run on cluster. Write job script for slurm to launch with
        # GNU parallel
        with cd(opt.output_dir):
            writeJobScript(opt, nRuns, cluster=True)

