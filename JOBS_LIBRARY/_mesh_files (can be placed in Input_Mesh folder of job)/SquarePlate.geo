h1 = 0.2;
h2 = 0.01;

L = 1;
H = 1;

Point(1) = {0  ,  -H/2, 0, h1};
Point(2) = {L  ,  -H/2, 0, h1};
Point(3) = {L  ,   H/2, 0, h1};
Point(4) = {0  ,   H/2, 0, h1};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};

Physical Surface(7) = {6};

Physical Line(8) = {1};
Physical Line(9) = {2};
Physical Line(10) = {3};
Physical Line(11) = {4};

Physical Point(12) = {1};
Physical Point(13) = {2};
Physical Point(14) = {3};
Physical Point(15) = {4};
