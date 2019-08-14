h1 = 0.1;

L = 1;
H = 1;

l = 0.5;
h = 0.001;

Point(1) = {0  , -H/2  , 0, h1};
Point(2) = {L  , -H/2  , 0, h1};
Point(3) = {L  ,  H/2  , 0, h1};
Point(4) = {0  ,  H/2  , 0, h1};

Point(5) = {0  , -h/2  , 0, h1};
Point(6) = {l  , -h/2  , 0, h1};
Point(7) = {l  ,  h/2  , 0, h1};
Point(8) = {0  ,  h/2  , 0, h1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 8};
Line(5) = {8, 7};
Line(6) = {7, 6};
Line(7) = {6, 5};
Line(8) = {5, 1};

Line Loop(9) = {3, 4, 5, 6, 7, 8, 1, 2};
Plane Surface(10) = {9};
Physical Surface(11) = {10};

Physical Line(12) = {1};
Physical Line(13) = {2};
Physical Line(14) = {3};
Physical Line(15) = {4};
Physical Line(16) = {8};

Physical Point(17) = {1};
Physical Point(18) = {2};
Physical Point(19) = {3};
Physical Point(20) = {4};
