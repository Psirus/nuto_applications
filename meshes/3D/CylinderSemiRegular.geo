Mesh.ElementOrder = 2;
Mesh.SecondOrderIncomplete = 1;

coarseness = 10;
Point(1) = {0, 0, 0, coarseness};
Point(2) = {31.0, 0, 0, coarseness};
Point(3) = {0, 31.0, 0, coarseness};
Point(4) = {0, -31.0, 0, coarseness};
Point(5) = {-31.0, 0, 0, coarseness};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 5};
Circle(3) = {5, 1, 4};
Circle(4) = {4, 1, 2};

Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};
Recombine Surface {6};

num[] = Extrude {0, 0, 186.0} { Surface{6}; Layers{31}; Recombine; };

Physical Volume (1) = {num[1]};
