Mesh.Algorithm = 6;
Mesh.HighOrderOptimize = 1;
Mesh.Optimize = 2;
Mesh.Smoothing = 2;
Mesh.SecondOrderIncomplete = 1;

h = 10;

r1 = 20.0;
Point(1) = {0, 0, 0, h};
Point(2) = {r1, 0, 0, h};
Point(3) = {0, r1, 0, h};
Point(4) = {0, -r1, 0, h};
Point(5) = {-r1, 0, 0, h};

r2 = 25.0;
Point(6) = {r2, 0, 0, h};
Point(7) = {0, r2, 0, h};
Point(8) = {0, -r2, 0, h};
Point(9) = {-r2, 0, 0, h};

r3 = 33.0;
Point(10) = {r3, 0, 0, h};
Point(11) = {0, r3, 0, h};
Point(12) = {0, -r3, 0, h};
Point(13) = {-r3, 0, 0, h};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 5};
Circle(3) = {5, 1, 4};
Circle(4) = {4, 1, 2};

Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 9};
Circle(7) = {9, 1, 8};
Circle(8) = {8, 1, 6};

Circle(9) = {10, 1, 11};
Circle(10) = {11, 1, 13};
Circle(11) = {13, 1, 12};
Circle(12) = {12, 1, 10};

Line Loop(13) = {2, 3, 4, 1};
Line Loop(14) = {6, 7, 8, 5};
Line Loop(15) = {10, 11, 12, 9};

Plane Surface(16) = {13};
Plane Surface(17) = {14, 13};
Plane Surface(18) = {15, 14};

Recombine Surface {16};
Recombine Surface {17};
Recombine Surface {18};

layers = 100.0/h;

concrete[] = Extrude {0, 0, 100} {
  Surface{16}; Layers{layers}; Recombine;
};

ceramics[] = Extrude {0, 0, 100} {
  Surface{17}; Layers{layers}; Recombine;
};

wool[] = Extrude {0, 0, 100} {
  Surface{18}; Layers{layers}; Recombine;
};

Physical Volume (1) = {concrete[1]};
Physical Volume (2) = {ceramics[1]};
Physical Volume (3) = {wool[1]};
