Mesh.Optimize = 1;
Mesh.Smoothing = 2;
Mesh.ElementOrder = 2;

xS = 0; xE = 100;
yS = 0; yE = 20;
meshSpecimen = 4;
// defines a box-shaped specimen 
// by start coordinates <xyz>S 
// and end coordinates  <xyz>E 

// points: 
p0 = newp; Point(p0) = {xS, yS, 0, meshSpecimen}; 
p1 = newp; Point(p1) = {xE, yS, 0, meshSpecimen}; 
p2 = newp; Point(p2) = {xE, yE, 0, meshSpecimen}; 
p3 = newp; Point(p3) = {xS, yE, 0, meshSpecimen}; 

// lines 
l0 = newreg; Line(l0) = {p0, p1}; 
l1 = newreg; Line(l1) = {p1, p2}; 
l2 = newreg; Line(l2) = {p2, p3}; 
l3 = newreg; Line(l3) = {p3, p0}; 

// lineloops and surfaces 
// the index says which coordinate is constant 
box = newreg; Line Loop(box) = { l0, l1, l2, l3}; 

volNr = newreg; 
Plane Surface(volNr) = {box};
Physical Surface(newreg) = volNr;
