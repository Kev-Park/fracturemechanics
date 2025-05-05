lc1 = 1;
lc2 = 0.2;
lc3 = 0.02;
h = 0.001;
a = 3;
r = 0.4;


Point(1) = {0,0,0,lc1};
Point(2) = {10,0,0,lc1};
Point(3) = {10,10,0,lc1};
Point(4) = {0,10,0,lc1};
Point(5) = {0,5+h,0,lc2};
Point(6) = {a,5,0,lc3};
Point(7) = {0,5-h,0,lc2};

Point(8) = {a+r*1.5,5,0,lc3};
Point(9) = {a+r,5+r,0,lc3};
Point(10) = {a+r,5-r,0,lc3};
Point(11) = {a,5+r,0,lc3};
Point(12) = {a,5-r,0,lc3};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,1};

Line Loop(1) = {1:7};
Plane Surface(1) = {1};
Point{8:12} In Surface{1};
