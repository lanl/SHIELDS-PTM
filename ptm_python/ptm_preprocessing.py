"""
This script file can be used to interpolate data from an SMWF run. It is assumed
that the data is in SWMF 3d idl format (binary OR ascii)
Binary is preferred due to the smaller size on disk
The ASCII tecplot format is also supported.
"""
import argparse
import glob
import os
import itertools
import numpy as np
# local package imports
try:
    from . import ptm_read
    from . import ptm_interpolate
except ImportError:
    # if running as a script in source directory
    # relative import won't work
    import ptm_read
    import ptm_interpolate

# Python 2 compatibility for checking string types
if not hasattr(__builtins__, "basestring"):
    basestring = (str, bytes)


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

    @classmethod
    def from_file(cls, filename):
        """Read PTM input format file
        """
        with open(filename, 'r') as fh:
            dims = fh.readline().strip().split()
            xvals = fh.readline().strip().split()
            yvals = fh.readline().strip().split()
            zvals = fh.readline().strip().split()
            datalines = fh.readlines()

        newobj = cls()
        # dimensions
        dims = [int(v) for v in dims]
        newobj.nx, newobj.ny, newobj.nz = dims
        newobj.x = np.array(xvals).astype(float)
        newobj.y = np.array(yvals).astype(float)
        newobj.z = np.array(zvals).astype(float)

        newobj.bx = np.empty((newobj.nx, newobj.ny, newobj.nz))
        newobj.by = np.empty((newobj.nx, newobj.ny, newobj.nz))
        newobj.bz = np.empty((newobj.nx, newobj.ny, newobj.nz))
        newobj.ex = np.empty((newobj.nx, newobj.ny, newobj.nz))
        newobj.ey = np.empty((newobj.nx, newobj.ny, newobj.nz))
        newobj.ez = np.empty((newobj.nx, newobj.ny, newobj.nz))

        #now populate (bx, by, bz) and (ex, ey, ez)
        for line in datalines:
            vals = line.strip().split()
            i, j, k = int(vals[0])-1, int(vals[1])-1, int(vals[2])-1
            b_e = [float(v) for v in vals[3:]]
            newobj.bx[i, j, k] = b_e[0]
            newobj.by[i, j, k] = b_e[1]
            newobj.bz[i, j, k] = b_e[2]
            newobj.ex[i, j, k] = b_e[3]
            newobj.ey[i, j, k] = b_e[4]
            newobj.ez[i, j, k] = b_e[5]

        return newobj


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


def tgrid_to_ascii(directory):
    """convert PTM tgrid file to ASCII format"""
    tgrid = np.fromfile(os.path.join(directory, 'tgrid.bin'))
    np.savetxt(os.path.join(directory, 'tgrid.dat'), tgrid)


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
    np.savetxt(os.path.join(dest_dir, 'tgrid.dat'), times)

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
        pf = PTMfields()
        pf.set_grid(xwant, ywant, zwant)
        pf.set_magnetic(res['Bx'], res['By'], res['Bz'])
        pf.set_electric(res['Ex'], res['Ey'], res['Ez'])
        pf.write_file(os.path.join(dest_dir, 'ptm_fields_{:04}.dat'.format(i1)))


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
