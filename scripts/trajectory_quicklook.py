import sys
import os
import argparse
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from ptm_python import ptm_tools

# This section only runs when launched as a command line script
# Either $> python scriptname [args]
# Or from ipython $> run -i scriptname [args]
if __name__ == '__main__':
    # Set up a basic argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--nparticle', dest='nparticle', type=int)
    parser.add_argument('-e', '--electrons', dest='electrons', action='store_true')
    parser.add_argument('fname')
    opt = parser.parse_args()

    # Load trajectory data
    data = ptm_tools.parse_trajectory_file(opt.fname)
    npart = 1 if opt.nparticle is None else opt.nparticle

    # Make a simple 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if npart < 1:
        # Plot trajectories for some random particle with access...
        keep = []
        for key in data:
            if not opt.electrons:
                size = 12
                if (np.linalg.norm(data[key][-1, 1:4]) >= 14.9) and (data[key][-1, 0] >= 180):
                    keep.append(key)
            else:
                if (np.linalg.norm(data[key][-1, 1:4]) >= 6.9) and (data[key][-1, 1] <= -2.5):
                    keep.append(key)
                size = np.min([45, len(keep)])
        usepart = np.random.choice(keep, size=size, replace=False)
        tts = []
        sit = data[1][0, 0]
        for pn in usepart:
            pldata = data[pn].copy()
            if not opt.electrons:
                ends = np.array([True] * len(data[pn])).astype(np.bool_)
                ends[:10] = False
                ends[-10:] = False
                pldata[ends] = np.nan
                ax.plot(pldata[:, 1], pldata[:, 2], pldata[:, 3])
            else:
                #ax.plot(pldata[-5:, 1], pldata[-5:, 2], pldata[-5:, 3])
                ax.plot(pldata[:-120, 1], pldata[:-120, 2], pldata[:-120, 3], lw=0.8, color=[0.5, 0.5, 0.5, 0.3])
                ax.plot(pldata[-120:, 1], pldata[-120:, 2], pldata[-120:, 3], alpha=0.6, lw=1.6)
                ax.plot([data[pn][-1, 1]], [data[pn][-1, 2]], [data[pn][-1, 3]],
                        linestyle='none', marker='o', ms=2, color='k')
                tts.append(sit - pldata[-1, 0])
        ax.plot([data[pn][0, 1]], [data[pn][0, 2]], [data[pn][0, 3]], linestyle='none', marker='o', color='k')
        if opt.electrons:
            print('Used {} particles for plot'.format(size))
            print('Min, max, and mean integration times of particles are:\n{}, {}, {}'.format(np.min(tts),
                                                                                              np.max(tts),
                                                                                              np.mean(tts)))
            ax.set_xlim([0, -7.2])
            ax.set_ylim([7.2, -7.2])
            ax.set_zlim([-3, 3])
    else:
        akima_x = scipy.interpolate.Akima1DInterpolator(1+data[npart][::-1, 0], data[npart][::-1, 1], axis=0)
        akima_y = scipy.interpolate.Akima1DInterpolator(1+data[npart][::-1, 0], data[npart][::-1, 2], axis=0)
        akima_z = scipy.interpolate.Akima1DInterpolator(1+data[npart][::-1, 0], data[npart][::-1, 3], axis=0)
        time = 1+np.linspace(data[npart][0, -1], data[npart][0, 0], len(data[npart])*5)
        px = akima_x(time)
        py = akima_y(time)
        pz = akima_z(time)
        ax.plot([data[npart][0, 1]], [data[npart][0, 2]], [data[npart][0, 3]], linestyle='none', marker='o')
        #ax.plot(data[npart][:, 1], data[npart][:, 2], data[npart][:, 3])
        ax.plot(px, py, pz)
    ax.set_xlabel('X$_{GSM}$ [R$_E$]')
    ax.set_ylabel('Y$_{GSM}$ [R$_E$]')
    ax.set_zlabel('Z$_{GSM}$ [R$_E$]')
    plt.show()
