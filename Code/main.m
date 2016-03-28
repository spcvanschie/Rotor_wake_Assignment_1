N = [150]; % number of annuli
lambda = [6,8,10]; % tip speed ratio [-]
U_inf = 10; % freestream velocity [m/s]
mu_min = 0.2; % spanwise start of blade [-]
rho = 1.225; % air density [kg/m^3]

% Optimisation parameters
pitch = 2; % blade pitch angle [deg]
twist = 14*(1-mu_local); % twist for each annulus [deg]
chordlength = 3*(1-mu_local)+1; % chord length for each annulus [m]

optimise = 1; % bogey statement to either use or ignore the optimisation part of the code

mu=cell(length(N),3);
a_all=cell(length(N),3);
a_tan_all=cell(length(N),3);
AoA_all=cell(length(N),3);
phi_all=cell(length(N),3);
C_t_all=cell(length(N),3);
C_n_all=cell(length(N),3);
Thrust_all = cell(length(N),3);
Torque_all = cell(length(N),3);
Cq_all = cell(length(N),3);
Power = ones(1,3);
for j = (1:length(N))
    for i = (1:3)
        % calculate interpolation splines for airfoil Cl and Cd
        [Clspline,Cdspline]=airfoil_liftdrag();

        % calculate geometrical parameters for each annulus
        [r,R,B,mu_min,mu_local,twist,chordlength,chordangle,omega,blade_solidity]= geometry(N(j),pitch,lambda(i),U_inf,mu_min,twist,chordlength);

        % calculate annulus characteristics, contains iteration loop for induction factors
        [W,phi,AoA,Cx,Cy,a_new,a_tan_new,Torque,C_torque,Thrust,Cp_all,P]=annulus_calc(rho,N(j),U_inf,r,R,omega,chordlength,chordangle,Clspline,Cdspline,blade_solidity,B,mu_local,lambda(i),mu_min,0,0);

    % writing all obtained variables to the memory    
    mu{j,i} = mu_local;
    a_all{j,i} = a_new;
    a_tan_all{j,i} = a_tan_new;
    AoA_all{j,i} = AoA;
    phi_all{j,i} = phi;
    C_t_all{j,i} = Glauert(a_new);
    C_n_all{j,i} = Cx;
    Thrust_all{j,i} = Thrust;
    Torque_all{j,i} = Torque;
    Cq_all{j,i} = C_torque;
    Power(1,i) = P;
    end
end

%% Blade optimisation section of code
% starting values of optimisation parameters
maxtwist = 14; % root twist angle [deg]
rootminustip = 3; % root chord length minus tip chord [m]
pitch = 2; % blade pitch angle [deg]

x0 = [maxtwist rootminustip pitch];
lambda = 8;

if optimise > 0
    C_t_design = 0.75;
    a = a_from_C_t(C_t_design)*ones(1,max(N));    
    
    twist = maxtwist*(1-mu_local); % twist for each annulus [deg]
    chordlength = rootminustip*(1-mu_local)+1; % chord length for each annulus [m]
    
    % calculate geometrical parameters for each annulus
    [r,R,B,mu_min,mu_local,twist,chordlength,chordangle,omega,blade_solidity]= geometry(N(j),pitch,lambda,U_inf,mu_min,twist,chordlength);

    % calculate annulus characteristics, contains iteration loop for induction factors
    [W,phi,AoA,Cx,Cy,a_new,a_tan_new,Torque,C_torque,Thrust,Cp,P]=annulus_calc(rho,N(j),U_inf,r,R,omega,chordlength,chordangle,Clspline,Cdspline,blade_solidity,B,mu_local,lambda,mu_min,optimise,a);

    
end

%% Plotting section of code
% figure axis ranges
axis_alpha = [0.2 1 -2 15];
axis_phi = [0.2 1 -2 30];
axis_a = [0.2 1 0 1];
axis_a_tan = [0.2 1 0 0.1];
axis_C_t = [0.2 1 0 1.2];
axis_C_n = [0.2 1 0 1.2];
axis_C_q = [0.2 1 0 1.5];

figure(1)
subplot(2,1,1)
plot(mu_local,AoA_all{length(N),1},mu_local,AoA_all{length(N),2},mu_local,AoA_all{length(N),3})
grid on
title('Angle of Attack')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
axis(axis_alpha)
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Angle of Attack \alpha [deg]')
subplot(2,1,2)
plot(mu_local,phi_all{length(N),1},mu_local,phi_all{length(N),2},mu_local,phi_all{length(N),3})
grid on
title('Inflow angle')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
axis(axis_phi)
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Inflow angle \phi [deg]')

figure(2)
subplot(2,1,1)
plot(mu_local,a_all{length(N),1},mu_local,a_all{length(N),2},mu_local,a_all{length(N),3})
grid on
title('Induction factor')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
axis(axis_a)
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Induction factor a [-]')
subplot(2,1,2)
plot(mu_local,a_tan_all{length(N),1},mu_local,a_tan_all{length(N),2},mu_local,a_tan_all{length(N),3})
grid on
title('Tangential induction factor')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
axis(axis_a_tan)
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Induction factor a'' [-]')

figure(3)
subplot(2,1,1)
plot(mu_local,C_t_all{length(N),1},mu_local,C_t_all{length(N),2},mu_local,C_t_all{length(N),3})
grid on
title('Thrust coefficient')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
axis(axis_C_t)
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Thrust coefficient C_{T} [-]')
subplot(2,1,2)
plot(mu_local,C_n_all{length(N),1},mu_local,C_n_all{length(N),2},mu_local,C_n_all{length(N),3})
grid on
title('Normal force coefficient')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
axis(axis_C_n)
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Normal force coefficient C_{n} [-]')

figure(4)
subplot(2,1,1)
plot(mu_local,Torque_all{length(N),1},mu_local,Torque_all{length(N),2},mu_local,Torque_all{length(N),3})
grid on
title('Torque')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
%axis(axis_C_t)
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Torque [N*m]')
subplot(2,1,2)
plot(mu_local,Cq_all{length(N),1},mu_local,Cq_all{length(N),2},mu_local,Cq_all{length(N),3})
grid on
title('Torque coefficient')
legend('\lambda = 6','\lambda = 8','\lambda = 10')
%axis(axis_C_q)
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Torque coefficient C_{q} [-]')

if length(N)>1
    figure(5)
    plot(mu{1,2},a_all{1,2},mu{2,2},a_all{2,2},mu{3,2},a_all{3,2})
    grid on
    title('Induction factor (\lambda = 8)')
    legend('N = 100','N = 300','N = 500')
    axis(axis_a)
    xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
    ylabel('Induction factor a [-]')
end