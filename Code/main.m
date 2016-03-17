N = 3; % number of annuli
pitch = 2; % blade pitch angle [deg]
lambda = 6; % tip speed ratio [-]
U_inf = 10; % freestream velocity [m/s]

% calculate interpolation splines for airfoil Cl and Cd
[Clspline,Cdspline]=airfoil_liftdrag();

% calculate geometrical parameters for each annulus
[r,R,B,mu_min,mu_local,twist,chordlength,chordangle,omega,blade_solidity]= geometry(N,pitch,lambda,U_inf);

% Initial guesses for a and a_tangential
a = 0.3*ones(1,N);
a_tan = 0.2*ones(1,N);

% calculate annulus characteristics, contains iteration loop for induction factors
[W,phi,AoA,Cx,Cy,a_new,a_tan_new]=annulus_calc(N,a,a_tan,U_inf,r,omega,chordangle,Clspline,Cdspline,blade_solidity);