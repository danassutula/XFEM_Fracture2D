
// 15946
// 62510
// 251546
// 1007777


// CASE #2
h0 = 0.10 * 1;
a = -6.1;
b = -3.6;

// CASE #3
// h0 = 0.075 * 8;
// a = -5.1;
// b = -3.4;

h1 = h0/16;
h2 = h0/16;

Point(1)  = {-10, -4, 0, 		h0};
Point(2)  = {-9, -4, 0, 		h0};
Point(3)  = {9, -4, 0, 			h0};
Point(4)  = {10, -4, 0, 		h0};
Point(5)  = {10, 4, 0, 			h0};
Point(6)  = {0, 4, 0, 			h0};
Point(7)  = {-10, 4, 0, 		h0};

Point(8)  = {-4, 2.75, 0, 		h1};
Point(9)  = {-4, 0.75, 0, 		h1};
Point(10) = {-4, -1.25, 0, 		h1};
Point(11) = {-4, -1.5, 0, 		h1};
Point(12) = {-3.75, -1.25, 0, 	h1};
Point(13) = {-4, -1.0, 0, 			h2};
Point(14) = {-4.25, -1.25, 0, 		h2};
Point(15) = {-4, 0.5, 0, 			h2};
Point(16) = {-4, 1, 0, 			h1};
Point(17) = {-3.75, 0.75, 0, 		h2};
Point(18) = {-4.25, 0.75, 0, 	h1};
Point(19) = {-4, 2.5, 0, 		h1};
Point(20) = {-4, 3, 0, 			h1};
Point(21) = {-3.75, 2.75, 0, 	h1};
Point(22) = {-4.25, 2.75, 0, 	h1};

Point(23) = {a, -4.0, 0, 	h1};
Point(24) = {b, -4.0, 0, 	h1};
Point(25) = {a,  4.0, 0, 	h1};
Point(26) = {b,  4.0, 0, 	h1};

Circle(8)  = {12, 10, 13};
Circle(9)  = {13, 10, 14};
Circle(10) = {14, 10, 11};
Circle(11) = {11, 10, 12};
Circle(12) = {17, 9, 16};
Circle(13) = {16, 9, 18};
Circle(14) = {18, 9, 15};
Circle(15) = {15, 9, 17};
Circle(16) = {21, 8, 20};
Circle(17) = {20, 8, 22};
Circle(18) = {22, 8, 19};
Circle(19) = {19, 8, 21};

Line(20) = {1, 2};
Line(21) = {2, 23};
Line(22) = {23, 24};
Line(23) = {24, 3};
Line(24) = {3, 4};
Line(25) = {4, 5};
Line(26) = {5, 6};
Line(27) = {6, 26};
Line(28) = {26, 25};
Line(29) = {25, 7};
Line(30) = {7, 1};
Line(31) = {23, 25};
Line(32) = {26, 24};

Line Loop(33) = {21, 31, 29, 30, 20};
Plane Surface(34) = {33};

Line Loop(35) = {22, -32, 28, -31};
Line Loop(36) = {18, 19, 16, 17};
Line Loop(37) = {13, 14, 15, 12};
Line Loop(38) = {9, 10, 11, 8};
Plane Surface(39) = {35, 36, 37, 38};

Line Loop(40) = {23, 24, 25, 26, 27, 32};
Plane Surface(41) = {40};

Characteristic Length {23} = h1;
Characteristic Length {24} = h1;
Characteristic Length {26} = h1;
Characteristic Length {25} = h1;

Physical Surface(42) = {34, 39, 41};

Physical Point(43) = {2};
Physical Point(44) = {3};
Physical Point(45) = {6};
