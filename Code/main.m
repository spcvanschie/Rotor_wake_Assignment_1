N = 3; % number of annuli
pitch = 2; % blade pitch angle [deg]
lambda = 6; % tip speed ratio [-]
U_inf = 10; % freestream velocity [m/s]

% Initial guesses for a and a_tangential
a = 0.3;
a_tan = 0.2;

% calculate geometrical parameters for each annulus
[R,B,mu_min,mu_local,twist,chordlength,omega]= geometry(N,pitch,lambda,U_inf);

% calculate interpolation splines for airfoil Cl and Cd
[Clspline,Cdspline]=airfoil_liftdrag();