/********************************************************************* 
 *
 *  Gmsh 
 * 
 *  circle inside a rectangular domain
 *
 *********************************************************************/

h1 = 2.00;
h2 = 0.50;

// rectangular plate:
xa = 0;     ya = 0; 
xb = 20;    yb = 10;

Point(1) = { xb, ya, 0 , h2};
Point(2) = { xb, yb, 0 , h2};
Point(3) = { xa, yb, 0 , h2};
Point(4) = { xa, ya, 0 , h2};

// inner circle 1
xc = 3; yc = 7; r = 2;

Point(5) = { xc,   yc,   0 , h2};
Point(6) = { xc+r, yc,   0 , h2};
Point(7) = { xc,   yc+r, 0 , h2};
Point(8) = { xc-r, yc,   0 , h2};
Point(9) = { xc,   yc-r, 0 , h2};

// inner circle 2
xc = 17; yc = 3; r = 2;

Point(10) = { xc,   yc,   0 , h2};
Point(11) = { xc+r, yc,   0 , h2};
Point(12) = { xc,   yc+r, 0 , h2};
Point(13) = { xc-r, yc,   0 , h2};
Point(14) = { xc,   yc-r, 0 , h2};

Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};
Circle(9) = {11, 10, 12};
Circle(10) = {12, 10, 13};
Circle(11) = {13, 10, 14};
Circle(12) = {14, 10, 11};

Line Loop(13) = {1, 2, 3, 4};
Line Loop(14) = {9, 10, 11, 12};
Line Loop(15) = {5, 6, 7, 8};
Plane Surface(16) = {13, 14, 15};


Physical Surface(17) = {16};

Physical Line(18) = {1};
Physical Line(19) = {2};
Physical Line(20) = {3};
Physical Line(21) = {4};
