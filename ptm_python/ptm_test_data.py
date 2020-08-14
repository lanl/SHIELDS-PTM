import os
import numpy as np
from scipy import interpolate, linalg
from matplotlib import pyplot as plt
import ptm_dipole
import ptm_preprocessing
import ptm_drift
import ptm_tools


def make_test_input(runid):
    """ """
    with open('ptm_input/ptm_parameters_{:04}.txt'.format(runid), 'w') as fh:
        fh.write('{:}\t\t\t\tRunID this needs to match the 4 digit tag\n'.format(runid))
        fh.write('{:}\t\t\t\tRand number seed\n'.format(2001))
        fh.write('{:}\t\t\t\tNumber of particles\n'.format(1))
        fh.write('{:}\t\t\t\tSpatial dimensions of the fields\n'.format(3))
        fh.write('{:}\t\t\t\tDirection of trace: -1=backwards,1=forwards\n'.format(1))
        fh.write('{:}\t\t\t\tIndex of the first data file\n'.format(1))
        fh.write('{:}\t\t\t\tIndex of the last data file\n'.format(2))
        fh.write('{:}\t\t\t\tTotal number of time files)\n'.format(2))
        fh.write('{:}\t\t\t\tCadence of input files (s)\n'.format(86400.0))
        fh.write('{:}\t\t\t\tCadence of output files (s)\n'.format(1.0))
        fh.write('{:}\t\t\t\tTime integrator, 1=RK4, 2=RKSuite\n'.format(1))
        fh.write('{:}\t\t\t\tEquations to solve: -1=drift,0=switch,1=orbit\n'.format(-1))
        fh.write('{:}\t\t\t\tGyrophase switching: 1=random,2=brute,3=gradient\n'.format(3))
        fh.write('{:}\t\t\t\tNumber of brute force points to search (angular res)\n'.format(36))
        fh.write('{:}\t\t\t\tCharge in multiples of fundamental\n'.format(1.0))
        fh.write('{:}\t\t\t\tMass in multiples of electron mass\n'.format(1836.0))
        fh.write('{:}\t\t\t\tLower limit of time integration (s)\n'.format(0.0))
        fh.write('{:}\t\t\t\tUpper limit of time integration (s)\n'.format(400.0))
        fh.write('{:}\t\t\t\tWrite trajectories in flux map mode\n'.format(0))
        fh.write('{:}\t\t\t\tInner boundary condition\n'.format(1))
        fh.write('{:}\t\t\t\tOuter boundary condition\n'.format(1))

    with open('ptm_input/dist_density_{:04}.txt'.format(runid), 'w') as fh:
        fh.write('{:}\t\t\t\tidens\n'.format(1))
        fh.write('{:}\t\t\t\tInitial X position (Re)\n'.format(5.0))
        fh.write('{:}\t\t\t\tInitial Y position (Re)\n'.format(0.0))
        fh.write('{:}\t\t\t\tInitial Z position (Re)\n'.format(0.0))

    with open('ptm_input/dist_velocity_{:04}.txt'.format(runid), 'w') as fh:
        fh.write('{:}\t\t\t\tidist\n'.format(1))
        fh.write('{:}\t\t\t\tE(keV)\n'.format(3000.0))
        fh.write('{:}\t\t\t\tPitch Angle (Degrees)\n'.format(90.0))
        fh.write('{:}\t\t\t\tPhase Angle (Degrees)\n'.format(180.0))


def make_test_data():
    """ """
    tgrid = np.array([0.0, 86400.0])
    xgrid = np.linspace(-8, 8,120)
    ygrid = np.linspace(-7, 7,110)
    zgrid = np.linspace(-6, 6,100)
    bdip = np.array([[[ptm_dipole.dipole_field([x, y, z]) for z in zgrid] for y in ygrid] for x in xgrid])

    ex = np.zeros([xgrid.size, ygrid.size, zgrid.size])
    ey = np.zeros_like(ex)
    ez = np.zeros_like(ex)

    pf = ptm_preprocessing.PTMfields()
    pf.set_grid(xgrid, ygrid, zgrid)
    pf.set_magnetic(bdip[:,:,:,0], bdip[:,:,:,1], bdip[:,:,:,2])
    pf.set_electric(ex, ey, ez)

    # Write the same data to two files, could replace the second call with shutil.copyfile
    pf.write_file('ptm_data/ptm_fields_0001.dat')
    pf.write_file('ptm_data/ptm_fields_0002.dat')

    # Create the time grid
    np.savetxt('ptm_data/tgrid.dat', tgrid)


def do_drift_comparison(runid, show=True):
    """
    """
    traj = ptm_tools.parse_trajectory_file('ptm_output/ptm_{:04}.dat'.format(runid))

    # Need to divide energy by 1000 b/c PTM natively uses keV
    td = ptm_drift.T_drift(traj[1][0,6]/1000,traj[1][0,7],linalg.norm(traj[1][0,1:4]),q=1,mc2=938)
    spl = interpolate.UnivariateSpline(traj[1][:,0],traj[1][:,1],k=4,s=0)
    rts = spl.derivative().roots()
    if len(rts) < 2:
        raise ValueError("Not enough data to determine period")
    elif len(rts) == 2:
        tp = 2*(rts[1]-rts[0])
    elif len(rts) > 2:
        tp = rts[2]-rts[0]

    # Plot x-position time series and highlight the expected and modeled drift periods
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 4))
    dlt = traj[1][:,1].max() - traj[1][:,1].min()
    ax.plot([td, td], [traj[1][:, 1].min() - 0.025*dlt,0], 'r:')
    ax.plot([tp, tp], [0,traj[1][:, 1].max() + 0.025*dlt], 'g:')
    ax.legend(['Hamlin', 'PTM'], loc='lower right')
    ax.plot(traj[1][:, 0], traj[1][:, 1], '-')
    ax.plot([0, traj[1][:, 0].max()], [0, 0], 'k', lw=0.5)
    ax.set_xlim([0, traj[1][:, 0].max()])
    ax.set_ylim([traj[1][:, 1].min() - 0.025*dlt, traj[1][:, 1].max() + 0.025*dlt])
    ax.set_title('Drift Periods: Analytical (Hamlin) = {:5.2f} s, Simulated (PTM) = {:5.2f} s'.format(td, tp))
    if show:
        plt.show()
    return td, tp


if __name__ == '__main__':
    import subprocess
    from ptm_tools import cd
    # check currdir and if called from inside SHIELDS-PTM change to basedir
    basename = 'SHIELDS-PTM'
    cwd = os.path.abspath(os.path.curdir)
    if basename in cwd and basename != os.path.split(cwd)[1]:
        basedir = os.path.join(cwd.split(basename), basename)
    # if not called from inside SHIELDS-PTM then exit
    elif basename not in cwd:
        raise Exception('{}: must be called from within {} directory tree, not {}'
                        .format(__file__, basename, cwd))
    else:
        basedir = cwd
    with cd(basedir):
        runid = 2
        print('Making dipole test input files...')
        make_test_input(runid)
        print('Completed input files')
        print('Making dipole test data files...')
        make_test_data()
        print('Completed data files')
        # Run PTM using runid
        result = subprocess.run(["./ptm", "2"], stdout=subprocess.PIPE)
        print(result.stdout.decode('utf-8'))
        print('\n\nCompleted PTM simulation')
    td, tp = do_drift_comparison(runid, show=False)
    # Test that estimated drift from PTM is within 0.5 seconds of analytic
    np.testing.assert_allclose(tp, td, atol=0.5)
