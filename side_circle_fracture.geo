lc1 = 1;
lc2 = 0.2;
lc3 = 0.02;
lc4 = 0.0002;
h = 0.001;
a = 3;
r = 0.4;


Point(1) = {0,0,0,lc1};
Point(2) = {10,0,0,lc1};
Point(3) = {10,10,0,lc1};
Point(4) = {0,10,0,lc1};
Point(5) = {0,5+h,0,lc2};
Point(6) = {a,5+h,0,lc4};
Point(7) = {a+h,5,0,lc4};
Point(8) = {a,5-h,0,lc4};
Point(9) = {0, 5-h, 0, lc2};
Point(10) = {a,5,0,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};

//Line(6) = {6,7};
//Line(7) = {7,8};

Circle(6) = {6, 10, 7};
Circle(7) = {7, 10, 8};
Line(8) = {8, 9};
Line(9) = {9, 1};

Line Loop(1) = {1:9};


Plane Surface(1) = {1};
Point{10} In Surface{1};
