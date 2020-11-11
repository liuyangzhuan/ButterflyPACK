rad=0.1;
ht=0.1;
cl1 = 1;
Point(1) = {rad, 0, 0, cl1};
Point(2) = {0, rad, 0, cl1};
Point(3) = {-rad, 0, 0, cl1};
Point(4) = {0, -rad, 0, cl1};
Point(5) = {0, 0, 0, cl1};
Point(6) = {rad, 0, ht, cl1};
Point(7) = {0, rad, ht, cl1};
Point(8) = {-rad, 0, ht, cl1};
Point(9) = {0, -rad, ht, cl1};
Point(10) = {0, 0, ht, cl1};

Circle(1) = {1, 5, 2};
Circle(2) = {2, 5, 3};
Circle(3) = {3, 5, 4};
Circle(4) = {4, 5, 1};
Circle(5) = {6, 10, 7};
Circle(6) = {7, 10, 8};
Circle(7) = {8, 10, 9};
Circle(8) = {9, 10, 6};
Line(9) = {1, 6};
Line(10) = {2, 7};
Line(11) = {3, 8};
Line(12) = {4, 9};

Line Loop(1) = {1,2,3,4};
Ruled Surface(1) = {1};
Line Loop(2) = {-5,-8,-7,-6};
Ruled Surface(2) = {2};
Line Loop(3) = {4,9,-8,-12};
Ruled Surface(3) = {3};
Line Loop(4) = {3,12,-7,-11};
Ruled Surface(4) = {4};
Line Loop(5) = {2,11,-6,-10};
Ruled Surface(5) = {5};
Line Loop(6) = {1,10,-5,-9};
Ruled Surface(6) = {6};


Physical Surface(1) = {1, 2, 3, 4, 5, 6};