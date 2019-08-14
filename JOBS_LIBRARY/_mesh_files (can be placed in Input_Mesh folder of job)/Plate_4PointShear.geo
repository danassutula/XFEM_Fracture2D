h = 5;
hf= 3;

Point(1) = {   0 ,   0 , 0, h};
Point(2) = {  20 ,   0 , 0, h};

Point(3) = { 240 ,   0 , 0, hf};
Point(4) = { 440 ,   0 , 0, h};

Point(5) = { 440 , 100 , 0, h};
Point(6) = { 420 , 100 , 0, h};

Point(7) = { 222.5 , 100 , 0, h};
Point(8) = { 222.5 ,  80 , 0, hf};
Point(9) = { 217.5 ,  80 , 0, hf};
Point(10)= { 217.5 , 100 , 0, h};

Point(11)= { 200 , 100 , 0, h};
Point(12)= {   0 , 100 , 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 1};

Line Loop(13) = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2};
Plane Surface(14) = {13};

Physical Surface(15) = {14};

Physical Point(16) = {2};
Physical Point(17) = {3};
Physical Point(18) = {6};
Physical Point(19) = {11};
