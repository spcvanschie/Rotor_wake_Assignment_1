N = 3; % number of annuli
pitch = 2; % blade pitch angle [deg]
lambda = 6; % tip speed ratio [-]
U_inf = 10; % freestream velocity [m/s]

% Initial guesses for a and a_tangential
a = 0.3*ones(1,N);
a_tan = 0.2*ones(1,N);

% calculate interpolation splines for airfoil Cl and Cd
[Clspline,Cdspline]=airfoil_liftdrag();

% calculate geometrical parameters for each annulus
[r,R,B,mu_min,mu_local,twist,chordlength,chordangle,omega]= geometry(N,pitch,lambda,U_inf);

% calculate annulus characteristics
[r,W,phi,AoA]=annulus_calc(a,a_tan,U_inf,mu_local,r,R,omega,chordangle);