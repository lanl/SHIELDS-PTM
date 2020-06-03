import sys
import os
import argparse
import numpy as np
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
    parser.add_argument('fname')
    opt = parser.parse_args()

    # Load trajectory data
    data = ptm_tools.parse_trajectory_file(opt.fname)
    npart = 1 if opt.nparticle is None else opt.nparticle

    # Make a simple 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(data[npart][:,1], data[npart][:,2], data[npart][:,3])
    ax.set_xlabel('X$_{GSM}$ [R$_E$]')
    ax.set_ylabel('Y$_{GSM}$ [R$_E$]')
    ax.set_zlabel('Z$_{GSM}$ [R$_E$]')
    plt.show()
