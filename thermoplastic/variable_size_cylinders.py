#!/usr/bin/env python3
import multiprocessing
import os
import subprocess
import numpy as np

Header = """
r = {0};
h = {1};
coarseness = {2};
"""

Template = """
Mesh.Algorithm3D = 2;
Mesh.Optimize = 4;
Mesh.ElementOrder = 2;

Point(1) = {0,0,0,coarseness};
Point(2) = {r,0,0,coarseness};
Point(3) = {0,r,0,coarseness};
Point(4) = {-r,0,0,coarseness};
Point(5) = {0,-r,0,coarseness};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};
Point{1} In Surface{6};

out[] = Extrude {0,0,h} {
  Surface{6};
};

Physical Volume(101) = {out[1]};
"""


def create_meshfile(r, h):
    coarseness = r / 10.0
    geofilename = 'Cylinder{0}.geo'.format(r)
    with open(geofilename, 'w') as geofile:
        geofile.write(Header.format(r, h, coarseness)+Template)
    meshfilename = 'Cylinder{0}.msh'.format(r)
    run = subprocess.Popen(['gmsh', '-3', '-order', '2', geofilename, '-o', meshfilename])
    run.wait()
    return meshfilename


def run_simulation(r):
    ratio = 3.0
    h = ratio * r
    meshfile = create_meshfile(r, h)
    resultdir = os.path.join('convergenceRadius', str(r))
    os.mkdir(resultdir)
    run = subprocess.Popen(['./HTLVariableSize', meshfile, resultdir, '500.0', str(r), str(h)])
    run.wait()


if __name__ == "__main__":
    radii = np.linspace(10, 60, 11)
    radii = np.append(radii, [90, 120, 150, 180, 210])
    pool = multiprocessing.Pool(16)
    pool.map(run_simulation, radii)
