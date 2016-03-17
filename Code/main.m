N = 50; % number of annuli
pitch = 2; % blade pitch angle [deg]
lambda = 10; % tip speed ratio [-]
U_inf = 10; % freestream velocity [m/s]

% calculate interpolation splines for airfoil Cl and Cd
[Clspline,Cdspline]=airfoil_liftdrag();

% calculate geometrical parameters for each annulus
[r,R,B,mu_min,mu_local,twist,chordlength,chordangle,omega,blade_solidity]= geometry(N,pitch,lambda,U_inf);

% calculate annulus characteristics, contains iteration loop for induction factors
[W,phi_all,AoA,Cx,Cy,a_new,a_tan_new]=annulus_calc(N,U_inf,r,omega,chordangle,Clspline,Cdspline,blade_solidity);

figure(1)
plot(mu_local,a_new)
grid on
title('Induction factor')
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Induction factor a [-]')