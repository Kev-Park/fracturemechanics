lc0 = 5;
lc1 = 1;
lc2 = 0.2;
lc3 = 0.02;
h = 0.001;
a = 3;
r = 0.4;


Point(1) = {0,0,0,lc0};
Point(2) = {10,0,0,lc0};
Point(3) = {10,10,0,lc0};
Point(4) = {0,10,0,lc0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1:4};

Point(5) = {5,5,0,lc0};
Point(6) = {4.5,5,0,lc0};
Point(7) = {5,4.5,0,lc0};
Point(8) = {5.5,5,0,lc0};
Point(9) = {5,5.5,0,lc0};

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};
Line Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1, 2};
Point{5} In Surface{1};
