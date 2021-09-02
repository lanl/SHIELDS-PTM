import argparse
import glob
import os

import numpy as np
import spacepy.datamodel as dm
import vtk

import ptm_python.ptm_preprocessing as pre


def get_fields_from_file(fname, inside_earth=1e-3):
    # Get Grid Field Data
    field = pre.PTMfields.from_file(fname)
    # Set field inside Earth to floor value, to avoid
    # very large values in output grid
    mxin = np.logical_and(field.x >= -1, field.x <= 1)
    myin = np.logical_and(field.y >= -1, field.y <= 1)
    mzin = np.logical_and(field.z >= -1, field.z <= 1)
    field.bx[mxin, myin, mzin] = inside_earth
    field.by[mxin, myin, mzin] = inside_earth
    field.bz[mxin, myin, mzin] = inside_earth
    return field

def make_vtk_grid(field):
    """Make a VTK UniformGrid object and set grid params"""
    grid = vtk.vtkUniformGrid()  # each block is a uniform grid
    grid_origin = [field.x.min(), field.y.min(), field.z.min()]
    dx = np.abs(field.x[2] - field.x[1])
    grid.SetOrigin(*grid_origin)
    grid.SetSpacing(dx, dx, dx)
    grid.SetDimensions(len(field.x),
                       len(field.y),
                       len(field.z))  # number of points in each direction
    return grid

def add_vtk_3vector(x1d, y1d, z1d, varname, grid):
    vararray = vtk.vtkDoubleArray()
    vararray.SetNumberOfComponents(3)
    vararray.SetNumberOfTuples(grid.GetNumberOfPoints())
    for ind in range(len(x1d)):
        vararray.SetTuple3(ind, x1d[ind], y1d[ind], z1d[ind])
    grid.GetPointData().AddArray(vararray)
    vararray.SetName(varname)

def write_grid_to_file(fname, grid, ascii=False):
    # set up writer for binary XML output
    writer = vtk.vtkXMLImageDataWriter()
    if ascii:
        writer.SetDataModeToAscii()
    writer.SetFileName(fname)
    writer.SetInputData(grid)
    writer.Write()


if __name__ == '__main__':
    # Set up a basic argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--magnetic', dest='magnetic',
                        action='store_true', default=False)
    parser.add_argument('-e', '--electric', dest='electric',
                        action='store_true', default=False)
    parser.add_argument('-a', '--ascii', dest='ascii',
                        action='store_true', default=False)
    parser.add_argument('fname')
    # TODO: handle wildcards in fname to process multiple files
    opt = parser.parse_args()
    if not (opt.magnetic or opt.electric):
        mywarn = "One or more of the fields options must be selected for writing."
        curropt = "Magnetic (-b): {}\nElectric (-e): {}".format(opt.magnetic, opt.electric)
        raise ValueError('\n'.join([mywarn, curropt]))
    outpath, infname = os.path.split(opt.fname)
    outfname = infname.split('.')[0]

    field = get_fields_from_file(opt.fname)
    grid = make_vtk_grid(field)
    if opt.magnetic:
        x1d = field.bx.flatten(order='F')
        y1d = field.by.flatten(order='F')
        z1d = field.bz.flatten(order='F')
        add_vtk_3vector(x1d, y1d, z1d, 'B', grid)
    if opt.electric:
        x1d = field.ex.flatten(order='F')
        y1d = field.ey.flatten(order='F')
        z1d = field.ez.flatten(order='F')
        add_vtk_3vector(x1d, y1d, z1d, 'E', grid)
    write_grid_to_file(os.path.join(outpath, '.'.join([outfname, 'vti'])),
                       grid, ascii=opt.ascii)
