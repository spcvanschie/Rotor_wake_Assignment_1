N = 400; % number of annuli
pitch = 2; % blade pitch angle [deg]
lambda = [6,8,10]; % tip speed ratio [-]
U_inf = 10; % freestream velocity [m/s]

a_all=cell(1,3);
a_tan_all=cell(1,3);
AoA_all=cell(1,3);
phi_all=cell(1,3);
a_all=cell(1,3);
for i = (1:3)
    % calculate interpolation splines for airfoil Cl and Cd
    [Clspline,Cdspline]=airfoil_liftdrag();

    % calculate geometrical parameters for each annulus
    [r,R,B,mu_min,mu_local,twist,chordlength,chordangle,omega,blade_solidity]= geometry(N,pitch,lambda(i),U_inf);

    % calculate annulus characteristics, contains iteration loop for induction factors
    [W,phi,AoA,Cx,Cy,a_new,a_tan_new]=annulus_calc(N,U_inf,r,omega,chordangle,Clspline,Cdspline,blade_solidity);
a_all{1,i} = a_new;
a_tan_all{1,i} = a_tan_new;
AoA_all{1,i} = AoA;
phi_all{1,i} = phi;
end

figure(1)
subplot(2,1,1)
plot(mu_local,a_all{1,1},mu_local,a_all{1,2},mu_local,a_all{1,3})
grid on
title('Induction factor')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Induction factor a [-]')
subplot(2,1,2)
plot(mu_local,a_tan_all{1,1},mu_local,a_tan_all{1,2},mu_local,a_tan_all{1,3})
grid on
title('Induction factor')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Induction factor a'' [-]')