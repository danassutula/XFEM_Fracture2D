hc = 2;
hf = 1;

r = 5/64;


Point(1) = {0, -2.5, 0, hc};
Point(2) = {5, -2.5, 0, hc};
Point(3) = {5,  2.5, 0, hc};
Point(4) = {0,  2.5, 0, hc};
Point(5) = {2,  0.0, 0, hf};
Point(6) = {3,  0.0, 0, hf};

Point(7) = {2-r,  0, 0, hf};
Point(8) = {  2, -r, 0, hf};
Point(9) = {2+r,  0, 0, hf};
Point(10)= {  2,  r, 0, hf};

Point(11)= {3-r,  0, 0, hf};
Point(12)= {  3, -r, 0, hf};
Point(13)= {3+r,  0, 0, hf};
Point(14)= {  3,  r, 0, hf};


Line(1) 	= {1, 2};
Line(2) 	= {2, 3};
Line(3)		= {3, 4};
Line(4) 	= {4, 1};

Circle(5) 	= {7, 5, 8};
Circle(6) 	= {8, 5, 9};
Circle(7) 	= {9, 5, 10};
Circle(8) 	= {10, 5, 7};
Circle(9) 	= {11, 6, 12};
Circle(10) 	= {12, 6, 13};
Circle(11) 	= {13, 6, 14};
Circle(12) 	= {14, 6, 11};


Line Loop(13) 	  = {3, 4, 1, 2};
Line Loop(14) 	  = {7, 8, 5, 6};
Line Loop(15) 	  = {12, 9, 10, 11};
Plane Surface(16) = {13, 14, 15};


Physical Point(101) = {1};
Physical Point(102) = {2};
Physical Point(103) = {3};
Physical Point(104) = {4};

Physical Line(201) = {1};
Physical Line(202) = {2};
Physical Line(203) = {3};
Physical Line(204) = {4};

Physical Surface(300) = {16};



