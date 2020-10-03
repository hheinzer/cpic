SetFactory("OpenCASCADE");

Sphere(1) = {0, 0, 0, 0.05, -Pi/2, Pi/2, 2*Pi};

Physical Surface(1) = {1, 2, 3};

Transfinite Curve {3, 1} = 10 Using Progression 1;
Transfinite Curve {6} = 20 Using Progression 1;
Transfinite Curve {5} = 10 Using Progression 1;
