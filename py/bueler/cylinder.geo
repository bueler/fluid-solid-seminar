// Define geometry of a cylinder (disc) in a rectangular box.
// Step 1.  Generate readable file cylinder.msh using Gmsh (https://gmsh.info/).
//   $ gmsh -2 cylinder.geo
// Step 2.  Run Navier-Stokes simulation; note cylinder.py reads cylinder.msh.
//   $ [activate Firedrake venv]
//   $ python3 cylinder.py

// box dimensions & target h
ly = 3.0;
minx = -3.0;
maxx = 12.0;
hbox = 0.5;

// cylinder dimensions & target h
rcyl = 1.0;
hcyl = 0.1;

// create geometry
Point(1) = {minx,  ly, 0, hbox};
Point(2) = {minx, -ly, 0, hbox};
Point(3) = {maxx, -ly, 0, hbox};
Point(4) = {maxx,  ly, 0, hbox};
Point(5) = { 0,  0, 0, hcyl};
Point(6) = { 1,  0, 0, hcyl};
Point(7) = {-1,  0, 0, hcyl};
Point(8) = { 0,  1, 0, hcyl};
Point(9) = { 0, -1, 0, hcyl};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Circle(5) = {8, 5, 6};
Circle(6) = {6, 5, 9};
Circle(7) = {9, 5, 7};
Circle(8) = {7, 5, 8};
Curve Loop( 9) = {1, 2, 3, 4};
Curve Loop(10) = {8, 5, 6, 7};
Plane Surface(1) = {9, 10};
Plane Surface(2) = {10};
Physical Curve("Upstream", 11) = {4};
Physical Curve("Downstream", 12) = {2};
Physical Curve("TopBottom", 13) = {1, 3};
Physical Curve("Circle", 14) = {8, 7, 6, 5};
Physical Surface("Domain", 3) = {1};
// Physical Surface("Disc", 4) = {2};
