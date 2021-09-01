import argparse
import numpy as np
from scipy import interpolate
import vtk
from vtk.util import numpy_support
from ptm_python import ptm_tools

if __name__ == '__main__':
    # Set up a basic argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--nparticle', dest='nparticle', type=int, default=-1)
    parser.add_argument('-r', '--random', dest='random', type=int, default=0)
    parser.add_argument('-s', '--smooth', dest='smooth', action='store_true')
    parser.add_argument('fname')
    opt = parser.parse_args()

    # Load trajectory data
    data = ptm_tools.parse_trajectory_file(opt.fname)

    if opt.random > 0:
        # Let's ignore named particle and pick some random particles that
        # have access
        keep = []
        for key in data:
            if (np.linalg.norm(data[key][-1, 1:4]) >= 14.9) and (data[key][-1, 0] >= 180):
                keep.append(key)
        size = np.min([opt.random, len(keep)])  # make sure we don't overflow
        usepart = np.random.choice(keep, size=size, replace=False)

    # For requested particle in the data file
    # TODO: need to loop over all particles and dump to combined file (XDMF?)
    if opt.nparticle >= 0 and opt.random <= 0:
        usepart = [opt.nparticle]

    for tnum in usepart:
        traj_pts = data[tnum]
        times = traj_pts[:, 0]  # time
        xyz_pts = traj_pts[:, 1:4]  # [x, y, z] position
        vperp = traj_pts[:, 4]  # field perpendicular velocity
        vpara = traj_pts[:, 5]  # field parallel velocity
        energy = traj_pts[:, 6]/1e3  # particle energy [MeV]
        alpha = traj_pts[:, 7]  # pitch angle
        if opt.smooth:
            akima_x = interpolate.Akima1DInterpolator(1+times[::-1], xyz_pts[::-1, 0], axis=0)
            akima_y = interpolate.Akima1DInterpolator(1+times[::-1], xyz_pts[::-1, 1], axis=0)
            akima_z = interpolate.Akima1DInterpolator(1+times[::-1], xyz_pts[::-1, 2], axis=0)
            akima_e = interpolate.Akima1DInterpolator(1+times[::-1], energy[::-1], axis=0)
            akima_a = interpolate.Akima1DInterpolator(1+times[::-1], alpha[::-1], axis=0)
            time = 1+np.linspace(times[-1], times[0], len(times)*5)
            px = akima_x(time)
            py = akima_y(time)
            pz = akima_z(time)
            pe = akima_e(time)
            pa = akima_a(time)
            xyz_pts = np.c_[px, py, pz]
        # Create a vtkPoints object and store the points in it
        points = vtk.vtkPoints()
        # Start by making a container for the cells...
        cells = vtk.vtkCellArray()
        # Create a polydata to store everything in
        polyData = vtk.vtkPolyData()
        vtk_energy = numpy_support.numpy_to_vtk(pe, deep=True, array_type=vtk.VTK_DOUBLE)
        vtk_energy.SetName('Energy')
        vtk_alpha = numpy_support.numpy_to_vtk(pe, deep=True, array_type=vtk.VTK_DOUBLE)
        vtk_alpha.SetName('Pitch Angle')
        for pt, en in zip(xyz_pts, pe):
            points.InsertNextPoint(pt)

        polyLine = vtk.vtkPolyLine()
        polyLine.GetPointIds().SetNumberOfIds(len(xyz_pts))
        for i in range(0, len(xyz_pts)):
            polyLine.GetPointIds().SetId(i, i)

        # Create a cell array to store the lines in and add the lines to it
        cells.InsertNextCell(polyLine)

        # Add the points to the dataset
        polyData.SetPoints(points)
        polyData.GetPointData().AddArray(vtk_energy)
        polyData.GetPointData().AddArray(vtk_alpha)
        polyData.GetPointData().SetActiveScalars('Energy')
        # Add the lines to the dataset
        polyData.SetLines(cells)

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetInputData(polyData)
        writer.SetDataModeToBinary()
        writer.SetFileName("ptm_particle{:04d}_E{:04d}.vtp".format(tnum, int(energy[0])))
        writer.Write()
