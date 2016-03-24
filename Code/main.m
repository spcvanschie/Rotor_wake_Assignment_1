N = 550; % number of annuli
pitch = 2; % blade pitch angle [deg]
lambda = [6,8,10]; % tip speed ratio [-]
U_inf = 10; % freestream velocity [m/s]
mu_min = 0.2; % spanwise start of blade [-]
rho = 1.225 % air density [kg/m^3]
delta_mu = (1-mu_min-(2*(1-mu_min)/N))/N; % width of each annulus [-]

a_all=cell(1,3);
a_tan_all=cell(1,3);
AoA_all=cell(1,3);
phi_all=cell(1,3);
C_t_all=cell(1,3);
C_n_all=cell(1,3);
Torque_all = cell(1,3);
Cq_all = cell(1,3);
for i = (1:3)
    % calculate interpolation splines for airfoil Cl and Cd
    [Clspline,Cdspline]=airfoil_liftdrag();

    % calculate geometrical parameters for each annulus
    [r,R,B,mu_min,mu_local,twist,chordlength,chordangle,omega,blade_solidity]= geometry(N,pitch,lambda(i),U_inf,mu_min);

    % calculate annulus characteristics, contains iteration loop for induction factors
    [W,phi,AoA,Cx,Cy,a_new,a_tan_new,Torque,C_torque]=annulus_calc(rho,N,U_inf,r,R,omega,chordlength,chordangle,Clspline,Cdspline,blade_solidity,B,mu_local,lambda(i),mu_min,delta_mu);
a_all{1,i} = a_new;
a_tan_all{1,i} = a_tan_new;
AoA_all{1,i} = AoA;
phi_all{1,i} = phi;
[C_t_all{1,i}] = Glauert(a_new);
C_n_all{1,i} = Cx;
Torque_all{1,i} = Torque;
Cq_all{1,i} = C_torque;
end

%% Plotting section of code
% figure axis ranges
axis_alpha = [0.2 1 0 15];
axis_phi = [0.2 1 0 30];
axis_a = [0.2 1 0 1];
axis_a_tan = [0.2 1 0 0.1];
axis_C_t = [0.2 1 0 1.2];
axis_C_n = [0.2 1 0 1.2];
axis_C_q = [0.2 1 0 1.5];

figure(1)
subplot(2,1,1)
plot(mu_local,AoA_all{1,1},mu_local,AoA_all{1,2},mu_local,AoA_all{1,3})
grid on
title('Angle of Attack')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
axis(axis_alpha)
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Angle of Attack \alpha [deg]')
subplot(2,1,2)
plot(mu_local,phi_all{1,1},mu_local,phi_all{1,2},mu_local,phi_all{1,3})
grid on
title('Inflow angle')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
axis(axis_phi)
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Inflow angle \phi [deg]')

figure(2)
subplot(2,1,1)
plot(mu_local,a_all{1,1},mu_local,a_all{1,2},mu_local,a_all{1,3})
grid on
title('Induction factor')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
axis(axis_a)
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Induction factor a [-]')
subplot(2,1,2)
plot(mu_local,a_tan_all{1,1},mu_local,a_tan_all{1,2},mu_local,a_tan_all{1,3})
grid on
title('Tangential induction factor')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
axis(axis_a_tan)
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Induction factor a'' [-]')

figure(3)
subplot(2,1,1)
plot(mu_local,C_t_all{1,1},mu_local,C_t_all{1,2},mu_local,C_t_all{1,3})
grid on
title('Thrust coefficient')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
axis(axis_C_t)
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Thrust coefficient C_{T} [-]')
subplot(2,1,2)
plot(mu_local,C_n_all{1,1},mu_local,C_n_all{1,2},mu_local,C_n_all{1,3})
grid on
title('Normal force coefficient')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
axis(axis_C_n)
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Normal force coefficient C_{n} [-]')

figure(4)
subplot(2,1,1)
plot(mu_local,Torque_all{1,1},mu_local,Torque_all{1,2},mu_local,Torque_all{1,3})
grid on
title('Torque')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
%axis(axis_C_t)
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Torque [N*m]')
subplot(2,1,2)
plot(mu_local,Cq_all{1,1},mu_local,Cq_all{1,2},mu_local,Cq_all{1,3})
grid on
title('Torque coefficient')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
%axis(axis_C_q)
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Torque coefficient C_{q} [-]')