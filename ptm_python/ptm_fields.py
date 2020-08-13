import os
import itertools
import numpy as np


class PTMfields(object):
    """Object for storing/writing PTM fields

    Example
    -------
    >>> pf = PTMfields()
    >>> pf.set_grid(xgrid, ygrid, zgrid)
    >>> pf.set_magnetic(bx, by, bz)
    >>> pf.set_electric(ex, ey, ez)
    >>> pf.write_file(os.path.join(directory, 'ptm_fields_{:04}.dat'.format(runid)))
    """
    def __init__(self):
        pass

    def set_grid(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.nx = np.size(x)
        self.ny = np.size(y)
        self.nz = np.size(z)

    def set_magnetic(self, bx, by, bz):
        bxs = np.shape(bx)
        bys = np.shape(by)
        bzs = np.shape(bz)

        if ~np.all([bxs == (self.nx, self.ny, self.nz),
                    bys == bxs, bzs == bxs]):
            raise ValueError('dimension mismatch in set_magnetic')

        self.bx = bx
        self.by = by
        self.bz = bz

    def set_electric(self, ex, ey, ez):
        exs = np.shape(ex)
        eys = np.shape(ey)
        ezs = np.shape(ez)

        if ~np.all([exs == (self.nx, self.ny, self.nz),
                    eys == exs, ezs == exs]):
            raise ValueError('dimension mismatch in set_magnetic')

        self.ex = ex
        self.ey = ey
        self.ez = ez

    def write_file(self, filename):
        """Write fields to PTM input format
        """
        with open(filename, 'w') as fh:
            fh.write('{:4} {:4} {:4}\n'.format(self.nx, self.ny, self.nz))
            fh.write((self.nx*'{:12.5e} ').format(*self.x).strip() + '\n')
            fh.write((self.ny*'{:12.5e} ').format(*self.y).strip() + '\n')
            fh.write((self.nz*'{:12.5e} ').format(*self.z).strip() + '\n')
            for i, j, k in itertools.product(range(self.nx),
                                             range(self.ny),
                                             range(self.nz)):
                fh.write(('{:4} {:4} {:4}'+(6*' {:12.5e}')).format(i+1, j+1, k+1,
                         self.bx[i, j, k],self.by[i, j, k],self.bz[i, j, k],
                         self.ex[i, j, k],self.ey[i, j, k],self.ez[i, j, k])+'\n')


def binary_to_xyz(directory, id):
    """Convert old-style binary PTM input files to new ASCII format

    Parameters
    ----------
    directory : string
        PTM_data directory with input files
    id : int
        Timestep number
    """
    xgrid = np.fromfile(os.path.join(directory, 'xgrid.bin'))
    ygrid = np.fromfile(os.path.join(directory, 'ygrid.bin'))
    zgrid = np.fromfile(os.path.join(directory, 'zgrid.bin'))
    nx, ny, nz = xgrid.size, ygrid.size, zgrid.size

    bx = np.fromfile(os.path.join(directory, 'bx3d_{:04}.bin'.format(id))).reshape([nx, ny, nz])
    by = np.fromfile(os.path.join(directory, 'by3d_{:04}.bin'.format(id))).reshape([nx, ny, nz])
    bz = np.fromfile(os.path.join(directory, 'bz3d_{:04}.bin'.format(id))).reshape([nx, ny, nz])
    ex = np.fromfile(os.path.join(directory, 'ex3d_{:04}.bin'.format(id))).reshape([nx, ny, nz])
    ey = np.fromfile(os.path.join(directory, 'ey3d_{:04}.bin'.format(id))).reshape([nx, ny, nz])
    ez = np.fromfile(os.path.join(directory, 'ez3d_{:04}.bin'.format(id))).reshape([nx, ny, nz])

    pf = PTMfields()
    pf.set_grid(xgrid, ygrid, zgrid)
    pf.set_magnetic(bx, by, bz)
    pf.set_electric(ex, ey, ez)
    pf.write_file(os.path.join(directory, 'ptm_fields_{:04}.dat'.format(id)))
