Mesh.Algorithm3D = 2;
Mesh.Optimize = 4;
Mesh.ElementOrder = 2;

coarseness = 31;
r = 31.0;
h = 186.0;

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

out[] = Extrude {0,0,h} {
  Surface{6};
};

Physical Volume(101) = {out[1]};
