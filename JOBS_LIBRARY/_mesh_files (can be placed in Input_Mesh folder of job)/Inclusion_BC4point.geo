/********************************************************************* 
 *
 *  Gmsh 
 * 
 *  circle inside a rectangular domain
 *
 *********************************************************************/

// defining some variables:
xa = 0; ya = 0; 
xb = 4; yb = 4; 
h1 = 0.75;

// POINTS
// outer domain
Point(1) = { xb, ya, 0 , h1};
Point(2) = { xb, yb, 0 , h1};
Point(3) = { xa, yb, 0 , h1};
Point(4) = { xa, ya, 0 , h1};

// inner circle 1
xc = 1; yc = 1; r = 0.5;

Point(5) = { xc,   yc,   0 , h1};
Point(6) = { xc+r, yc,   0 , h1};
Point(7) = { xc,   yc+r, 0 , h1};
Point(8) = { xc-r, yc,   0 , h1};
Point(9) = { xc,   yc-r, 0 , h1};

// inner circle 2
xc = 2.5; yc = 2.5; r = 0.8;

Point(10) = { xc,   yc,   0 , h1};
Point(11) = { xc+r, yc,   0 , h1};
Point(12) = { xc,   yc+r, 0 , h1};
Point(13) = { xc-r, yc,   0 , h1};
Point(14) = { xc,   yc-r, 0 , h1};

// LINES AND ARCS
// inner circle 1
Circle(1) = {6,5,7};
Circle(2) = {7,5,8};
Circle(3) = {8,5,9};
Circle(4) = {9,5,6};

// inner circle 2
Circle(5) = {11,10,12};
Circle(6) = {12,10,13};
Circle(7) = {13,10,14};
Circle(8) = {14,10,11};

// exterior
Line(9)  = {1,2};
Line(10) = {2,3};
Line(11) = {3,4};
Line(12) = {4,1};

//SURFACES
Line Loop(9)  = {1,2,3,4};     // interior
Line Loop(10) = {5,6,7,8};     // interior
Line Loop(11) = {9,10,11,12};  // exterior

Plane Surface(21) = {9};       // interior
Plane Surface(22) = {10};      // interior
Plane Surface(23) = {11,9,10}; // exterior, interior

Color Blue  { Surface{ 21:22 }; }
Color Red   { Surface{ 23 }; }

Physical Point(40) = {4};
Physical Point(41) = {1};
Physical Point(42) = {2};
Physical Point(43) = {3};
