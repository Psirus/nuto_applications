Mesh.Algorithm = 6;
Mesh.HighOrderOptimize = 1;
Mesh.Optimize = 2;
Mesh.Smoothing = 2;
Function MySphere 
  meshCircleR = meshCircle; 
  p1 = newp; Point(p1) = {xC,  yC,  zC,  meshCircleR} ; 
  p2 = newp; Point(p2) = {xC+R,yC,  zC,  meshCircleR} ; 
  p3 = newp; Point(p3) = {xC,  yC+R,zC,  meshCircleR} ; 
  p4 = newp; Point(p4) = {xC,  yC,  zC+R,meshCircleR} ; 
  p5 = newp; Point(p5) = {xC-R,yC,  zC,  meshCircleR} ; 
  p6 = newp; Point(p6) = {xC,  yC-R,zC,  meshCircleR} ; 
  p7 = newp; Point(p7) = {xC,  yC,  zC-R,meshCircleR} ; 
 
  c1  = newreg; Circle(c1)  = {p2,p1,p7}; c2  = newreg; Circle(c2)  = {p7,p1,p5}; 
  c3  = newreg; Circle(c3)  = {p5,p1,p4}; c4  = newreg; Circle(c4)  = {p4,p1,p2}; 
  c5  = newreg; Circle(c5)  = {p2,p1,p3}; c6  = newreg; Circle(c6)  = {p3,p1,p5}; 
  c7  = newreg; Circle(c7)  = {p5,p1,p6}; c8  = newreg; Circle(c8)  = {p6,p1,p2}; 
  c9  = newreg; Circle(c9)  = {p7,p1,p3}; c10 = newreg; Circle(c10) = {p3,p1,p4}; 
  c11 = newreg; Circle(c11) = {p4,p1,p6}; c12 = newreg; Circle(c12) = {p6,p1,p7}; 
 
  l1 = newreg; Line Loop(l1) = {c5,c10,c4};   Ruled Surface(newreg) = {l1}; 
  l2 = newreg; Line Loop(l2) = {c9,-c5,c1};   Ruled Surface(newreg) = {l2}; 
  l3 = newreg; Line Loop(l3) = {c12,-c8,-c1}; Ruled Surface(newreg) = {l3}; 
  l4 = newreg; Line Loop(l4) = {c8,-c4,c11};  Ruled Surface(newreg) = {l4}; 
  l5 = newreg; Line Loop(l5) = {-c10,c6,c3};  Ruled Surface(newreg) = {l5}; 
  l6 = newreg; Line Loop(l6) = {-c11,-c3,c7}; Ruled Surface(newreg) = {l6}; 
  l7 = newreg; Line Loop(l7) = {-c2,-c7,-c12};Ruled Surface(newreg) = {l7}; 
  l8 = newreg; Line Loop(l8) = {-c6,-c9,c2};  Ruled Surface(newreg) = {l8}; 
   
  theLoops[t] = newreg; 
 
  Surface Loop(theLoops[t]) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1}; 
 
  thehole = newreg; 
  Volume(thehole) = theLoops[t]; 
  theAggregates[t] = thehole; 
 
Return 
 
radius = 31;
height = 93;
meshSpecimen = 5;
// defines a cylinder-shaped specimen

// base cirle
p0 = newp; Point(p0) = {0, 0, 0, meshSpecimen};
p1 = newp; Point(p1) = {0, radius, 0, meshSpecimen};
p2 = newp; Point(p2) = {0, -radius, 0, meshSpecimen};
c0 = newreg; Circle(c0) = {p1, p0, p2};
c1 = newreg; Circle(c1) = {p2, p0, p1};
baseLoop = newreg; Line Loop(baseLoop) = {c0, c1};
baseSurface = newreg; Plane Surface(baseSurface) = {baseLoop};

// top cirle
p3 = newp; Point(p3) = {0, 0, height, meshSpecimen};
p4 = newp; Point(p4) = {0, radius, height, meshSpecimen};
p5 = newp; Point(p5) = {0, -radius, height, meshSpecimen};
c2 = newreg; Circle(c2) = {p4, p3, p5};
c3 = newreg; Circle(c3) = {p5, p3, p4};
topLoop = newreg; Line Loop(topLoop) = {c2, c3};
topSurface = newreg; Plane Surface(topSurface) = {topLoop};

// barrell surfaces
verticalOne = newreg; Line(verticalOne) = {p2, p5};
verticalTwo = newreg; Line(verticalTwo) = {p1, p4};
barrelLoopOne = newreg; Line Loop(barrelLoopOne) = {c0, verticalOne, -c2, -verticalTwo};
barrelLoopTwo = newreg; Line Loop(barrelLoopTwo) = {c1, verticalTwo, -c3, -verticalOne};
barrelSurfaceOne = newreg; Ruled Surface(barrelSurfaceOne) = {barrelLoopOne};
barrelSurfaceTwo = newreg; Ruled Surface(barrelSurfaceTwo) = {barrelLoopTwo};
// surface loop
theLoops[0] = newreg;
Surface Loop(theLoops[0]) = {barrelSurfaceOne, barrelSurfaceTwo, baseSurface, topSurface};
theBox = newreg;
Volume(theBox) = theLoops[0];


meshCircle = 5; 
t = 1;
xC = -8.23641; yC = -16.4615; zC = 80.407;
R = 7.77757; 
Call MySphere; 
 
 
t = 2;
xC = -16.1444; yC = 9.17825; zC = 41.0489;
R = 7.61355; 
Call MySphere; 
 
 
t = 3;
xC = -9.65526; yC = 16.1632; zC = 12.1726;
R = 7.35724; 
Call MySphere; 
 
 
t = 4;
xC = -18.2605; yC = -4.87431; zC = 61.0423;
R = 7.28452; 
Call MySphere; 
 
 
t = 5;
xC = 2.54669; yC = -18.8407; zC = 11.9884;
R = 7.1727; 
Call MySphere; 
 
 
t = 6;
xC = -17.9323; yC = -6.47671; zC = 11.9482;
R = 7.11869; 
Call MySphere; 
 
 
t = 7;
xC = -12.082; yC = -14.7913; zC = 40.0768;
R = 7.08603; 
Call MySphere; 
 
 
t = 8;
xC = 5.16277; yC = 18.832; zC = 45.9619;
R = 6.6578; 
Call MySphere; 
 
 
t = 9;
xC = 18.3702; yC = 6.83399; zC = 81.5994;
R = 6.58425; 
Call MySphere; 
 
 
t = 10;
xC = -0.273534; yC = -0.39187; zC = 26.0781;
R = 6.58023; 
Call MySphere; 
 
 
t = 11;
xC = 11.0816; yC = -16.2703; zC = 39.612;
R = 6.49843; 
Call MySphere; 
 
 
t = 12;
xC = -10.1847; yC = 16.9533; zC = 62.576;
R = 6.40721; 
Call MySphere; 
 
 
t = 13;
xC = -18.0031; yC = 6.70148; zC = 81.6857;
R = 6.35492; 
Call MySphere; 
 
 
t = 14;
xC = 16.793; yC = 2.55858; zC = 59.4163;
R = 6.33099; 
Call MySphere; 
 
 
t = 15;
xC = 0.227316; yC = 19.9863; zC = 81.9871;
R = 6.19719; 
Call MySphere; 
 
 
t = 16;
xC = 14.3795; yC = -13.9271; zC = 74.2517;
R = 6.16637; 
Call MySphere; 
 
 
t = 17;
xC = 0.161419; yC = -20.0257; zC = 58.6752;
R = 6.157; 
Call MySphere; 
 
 
t = 18;
xC = 19.6046; yC = 4.15201; zC = 37.0628;
R = 6.14436; 
Call MySphere; 
 
 
t = 19;
xC = 19.6299; yC = -4.57346; zC = 17.0795;
R = 6.02901; 
Call MySphere; 
 
 
t = 20;
xC = 13.2838; yC = 15.182; zC = 10.8272;
R = 6.01169; 
Call MySphere; 
 
 
volNr = newreg; 
Volume(volNr) = {theLoops[]};
Physical Volume(newreg) = volNr;
Physical Volume(newreg) = {theAggregates[]};
