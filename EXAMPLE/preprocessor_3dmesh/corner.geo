sph_rad=1;
cl1 = 1;
Point(1) = {sph_rad, 0, 0, cl1};
Point(2) = {0, sph_rad, 0, cl1};
Point(3) = {0, 0, sph_rad, cl1};
Point(4) = {0, 0, 0, cl1};

Line(1) = {4, 1};
Line(2) = {4, 2};
Line(3) = {4, 3};
Line(4) = {1, 3};
Line(5) = {3, 2};
Line(6) = {2, 1};


Line Loop(1) = {4,-3,1};
Ruled Surface(1) = {1};
Line Loop(2) = {3,5,-2};
Ruled Surface(2) = {2};
Line Loop(3) = {-1,2,6};
Ruled Surface(3) = {3};

Physical Surface(1) = {1, 2, 3};
