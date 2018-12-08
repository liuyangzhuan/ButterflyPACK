sph_rad=1;
cl1 = 1;
Point(1) = {sph_rad, 0, 0, cl1};
Point(2) = {0, sph_rad, 0, cl1};
Point(3) = {-sph_rad, 0, 0, cl1};
Point(4) = {0, -sph_rad, 0, cl1};
Point(5) = {0, 0, sph_rad, cl1};
Point(7) = {0, 0, 0, cl1};
Circle(1) = {1, 7, 2};
Circle(2) = {2, 7, 3};
Circle(3) = {3, 7, 4};
Circle(4) = {4, 7, 1};
Circle(5) = {1, 7, 5};
Circle(6) = {2, 7, 5};
Circle(7) = {3, 7, 5};
Circle(8) = {4, 7, 5};

Line Loop(1) = {1,6,-5};
Ruled Surface(1) = {1};
Line Loop(2) = {2,7,-6};
Ruled Surface(2) = {2};
Line Loop(3) = {3,8,-7};
Ruled Surface(3) = {3};
Line Loop(4) = {4,5,-8};
Ruled Surface(4) = {4};

Physical Surface(1) = {1, 2, 3, 4};

