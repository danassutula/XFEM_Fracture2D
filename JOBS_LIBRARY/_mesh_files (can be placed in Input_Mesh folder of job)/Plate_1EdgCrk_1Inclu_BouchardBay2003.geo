h1 = 0.2;
h2 = 0.2;

r = 1/3;


Point(1) = {0  , -1  , 0, h1};
Point(2) = {1  , -1  , 0, h1};
Point(3) = {1  ,  1  , 0, h1};
Point(4) = {0  ,  1  , 0, h1};

Point(5) = {0.5  , -0.5  , 0, h2};
Point(6) = {0.5-r, -0.5  , 0, h2};
Point(7) = {0.5  , -0.5-r, 0, h2};
Point(8) = {0.5+r, -0.5  , 0, h2};
Point(9) = {0.5  , -0.5+r, 0, h2};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};


Line Loop(9) = {2, 3, 4, 1};
Line Loop(10) = {7, 8, 5, 6};

Plane Surface(11) = {9, 10};
Plane Surface(12) = {10};


Physical Surface(13) = {11};
Physical Surface(14) = {12};

Physical Point(15) = {1};
Physical Point(16) = {2};
Physical Point(17) = {3};
Physical Point(18) = {4};

Physical Line(19) = {1};
Physical Line(20) = {2};
Physical Line(21) = {3};
Physical Line(22) = {4};

