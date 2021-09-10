#!/usr/bin/env python
"""
SHIELDS-PTM run setup script for SEP access
"""
import argparse
import glob
import os
import re
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
    indir = os.path.join(opt.input_dir, 'ptm_data')
    # get timesteps available in ptm_data directory
    gterm = os.path.join(indir, 'ptm_fields*dat')
    cands = glob.glob(gterm)
    cands = [re.search('ptm_fields_(\d{4})', cc) for cc in cands]
    cands = [int(cc.group().split("_")[-1]) for cc in cands]
    step1 = min(cands)
    stepN = max(cands)
    tgrid = np.loadtxt(os.path.join(indir, 'tgrid.dat'))
    # break energy & alpha space over number of runs...
    emin = 5.0*1e3  # 5MeV lower limit
    emax = 800.0*1e3  # 800MeV upper limit
    nenergy = 300
    pamin = 0.1
    pamax = 179.9
    nalpha = 67
    energy_arr = np.logspace(np.log10(emin), np.log10(emax), num=nenergy)
    per_run = nenergy//nruns
    # pitch angles evenly spaced around unit sphere
    # TODO: with one phase angle per (E,a) pair, coverage is relatively sparse
    #       perp. to the magnetic field. Can correct with, e.g., cosine scaling
    #       of alpha, but could switch to a method for evenly distributing points
    #       on the sphere (e.g., https://gist.github.com/dinob0t/9597525)
    #       This would then necessitate a different way of postprocessing the
    #       fluxmap
    alpha_arr = np.linspace(pamin, pamax, num=nalpha)
    # now reorder energy arrays randomly
    np.random.shuffle(energy_arr)

    for runno in range(nruns):
        # for each run, use a subset of the energies
        imin = per_run*runno
        imax = per_run*(runno+1)  # exclusive for use in ranges
        imin_a = 0  # will need to change for splits across alpha
        imax_a = nalpha
        earr = {'n_e': per_run if runno != nruns-1 else len(energy_arr[imin:]),
                'n_a': nalpha,
                }
        # build energy and pitch angle arrays and write to file
        # TODO: to use different pitch angles the different alpha sets
        #       would need to be paired with copies of the energy sets
        #       E.g. Energies are [1, 3, 5] and [2, 4]
        #            alphas are [0, 45, 90] and [22.5, 67.5]
        #       So E=[1, 3, 5] would need to be run with both a=[0, 45, 90] and a=[22.5, 67.5]
        runenergies = np.array(energy_arr[imin:imax])
        runalphas = np.array(alpha_arr[imin_a:imax_a])
        runenergies.tofile(os.path.join(indir, 'energies_{:04d}.bin'.format(runid+runno)))
        runalphas.tofile(os.path.join(indir, 'pitchangles_{:04d}.bin'.format(runid+runno)))

        # start input file setup
        setup = ptm_input.ptm_input_creator(runid=runid+runno,
                                            idensity=1,  # Point source
                                            )

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
                             ifirst=step1,  # number of first time
                             ilast=stepN,  # number of last time
                             ntot=len(tgrid),  # number of steps
                             charge=1.0, mass=1837.0,  # protons
                             dtin=opt.timestep,  # set delta between fields snapshots
                             dtout=0.1,  # 0.1s output written
                             thi=opt.starttime,  # max time ("start" for backwards trace)
                             tlo=opt.starttime-300,  # min time ("end" for backwards trace)
                             itraj=1,  # 1=write trajectories when in fluxmap mode
                             ibound=3,  # 3=radial distance stopping criterion
                             rbound=15.0,  # distance to stopping surface
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
    # break energy space over number of runs...
    emin = 100.0*1e-3  # 100eV lower limit
    emax = 100.0  # 200eV upper limit
    nenergy = 90
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
                             pamin=0,
                             pamax=90,
                             xsource=-12.0
                             )

        # Set tracing params
        setup.set_parameters(iswitch=0,  # switch between GC and full orbit
                             itrace=-1,  # backwards in time
                             charge=-1.0, mass=1.0,  # electrons
                             dtin=opt.timestep,  # set delta between fields snapshots
                             dtout=1.0, # 1.0s output written
                             thi=600.0, # max time of 600s
                             itraj=1,  # 1=write trajectories when in fluxmap mode
                             ibound=3,  # 2=plane at fixed X; 3=radial bound
                             rbound=-6.6,  # distance to stopping surface
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
        nRuns = 10
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

