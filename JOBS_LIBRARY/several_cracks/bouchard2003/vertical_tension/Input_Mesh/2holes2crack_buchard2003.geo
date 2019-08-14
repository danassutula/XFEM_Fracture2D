
hc = 0.2;
hf = 0.2;

r = 2;


Point(1) = {-10, -5, 0, hc};
Point(2) = { 10, -5, 0, hc};
Point(3) = { 10,  5, 0, hc};
Point(4) = {-10,  5, 0, hc};
Point(5) = { -7,  2, 0, hf};
Point(6) = {  7, -2, 0, hf};

Point(7) = {-7-r,  2, 0, hf};
Point(8) = {  -7,2+r, 0, hf};
Point(9) = {-7+r,  2, 0, hf};
Point(10)= {  -7,2-r, 0, hf};

Point(11)= {7-r,  -2, 0, hf};
Point(12)= {  7,-2-r, 0, hf};
Point(13)= {7+r,  -2, 0, hf};
Point(14)= {  7,-2+r, 0, hf};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Circle(5) = {10, 5, 9};
Circle(6) = {9, 5, 8};
Circle(7) = {8, 5, 7};
Circle(8) = {7, 5, 10};

Circle(9) = {12, 6, 13};
Circle(10) = {13, 6, 14};
Circle(11) = {14, 6, 11};
Circle(12) = {11, 6, 12};

Line Loop(13) = {10, 11, 12, 9};
Line Loop(14) = {6, 7, 8, 5};
Line Loop(15) = {3, 4, 1, 2};

Plane Surface(16) = {13, 14, 15};
