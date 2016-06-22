clear variables;
close all;
%% Rotor & flow parameter declaration
[Clspline,Cdspline,alphadata,Cldata,Cddata]=airfoil_liftdrag();
N = 3; % number of rotor blades
rho = 1.225; % air density [kg/m3]
span = 50; % span of the rotor blade [m]
blade_root = 0.2*span; % location of blade root [m]
rootchord = 4; % root chord length (at mu = 0) [m]
tipchord = 2; % tip chord length [m]
bladepitch = -1.5; % blade pitch angle [deg]
roottwist = 1; % root twist angle [deg]
tiptwist = 0; % tip twist angle [deg]
u_inf = [10;0;0]; % freestream velocity vector [m/s]
tsr = 8; % Tip Speed Ratio [-]
omega = tsr*u_inf(1)/span; % rotor rotational rate in rad/s

%% Discretisation parameter declaration
n = 160; % number of collocation points [-]
dt = 0.01; % time step [s]
max_timestep = 200; % number of time steps [-]
timestep = 0; % initial time step number [-]
rotation_angle = 0; % initial azimuthal blade rotation angle [rad], defined as 0 along the y-axis
eps = 10^-6; % vortex core epsilon [m]

%% Calculation of initial collocation point coordinates
cp_init_y = cosspace(blade_root,span,n+2,0);
% Removal of root and tip points from CP y-locations
cp_init_y = cp_init_y(cp_init_y~=blade_root);
cp_init_y = cp_init_y(cp_init_y~=span);

cp_init_z = zeros(1,n);

%% Generation of initial vortex element endpoint coordinates
vort_end_y = endpoints(blade_root,span,cp_init_y);
vort_end_z = zeros(1,length(vort_end_y));

vort_end_spanwise = vort_end_y; % radial position of every vortex element endpoint

%% Calculating (constant) geometrical parameters for every discretised element
chord = chordlength(span,rootchord,tipchord,cp_init_y);
alpha = blade_angle(span,roottwist,tiptwist,bladepitch,cp_init_y);
surface_area = surfarea(chord,vort_end_y);
panel_width = surface_area./chord;
cp_spanwise = cp_init_y; % radial position of every collocation point
cp_mu = cp_spanwise./span;

%% Adjustment of collocation point x-coordinates and trailing edge vortex corners
cp_init_x = 0.75.*chord;
cp_init = [cp_init_x; cp_init_y; cp_init_z];

vort_end_x = 0.25*chordlength(span,rootchord,tipchord,vort_end_y); % the vortex elements are located in the yz-plane
vort_end_init = [vort_end_x; vort_end_y; vort_end_z];

TE_x = chordlength(span,rootchord,tipchord,vort_end_y);
TE = [TE_x; vort_end_y; vort_end_z];

%% Calculation of initial collocation point normal vectors
normalvectors = normalvect(alpha,rotation_angle);

u_inf_repmat = repmat(u_inf,1,length(normalvectors(1,:)));

%% Calculating initial influence coefficient matrix & freestream and blade rotation influence vector
A_ind_init = induction_factors(normalvectors,cp_init,vort_end_init,TE,eps,0);
u_rotor = rotationvect(omega,rotation_angle,cp_spanwise,normalvectors);
Q_inf_init = freestreamvect(normalvectors,u_inf_repmat,u_rotor);
%Q_rot = rotationvect(omega,rotation_angle,cp_spanwise,normalvectors);

%% Solving the initial linear system of equations for the vortex circulations
circulations = -A_ind_init\(Q_inf_init);

%% Time marching section for t > 0

Thrust_convergence_ll = zeros(max_timestep,1);
Q_ind_start = zeros(length(chord),max_timestep-1);
%% Start of time marching
for timestep = 1:max_timestep
    %% Iteration loop
    
    cp_coords_blade = bladerotation(rotation_angle,cp_spanwise,cp_init);
    vort_end_blade = bladerotation(rotation_angle,vort_end_spanwise,vort_end_init);
    TE_blade = bladerotation(rotation_angle,vort_end_spanwise,TE);
    normalvectors_upd = normalvect(alpha,rotation_angle);
    u_rotor_upd = rotationvect(omega,rotation_angle,cp_spanwise,normalvectors_upd);
    
    
    A_ind_upd = induction_factors(normalvectors_upd,cp_coords_blade,vort_end_blade,TE_blade,eps,0);
    [Q_inf,u_combined] = freestreamvect(normalvectors_upd,u_inf_repmat,u_rotor_upd);
    
    if timestep < 2
        Q_ind = zeros(length(cp_spanwise),1);
    else
        [Q_ind, Q_coefs] = wakeinduction(Q_ind_start,cp_coords_blade,normalvectors_upd,TE_convected,TE_blade,circulation_history,eps,0);
    end
    
    circulation_blade = -A_ind_upd\(Q_inf + Q_ind);
%     for i = 1:length(circulation_blade);
%         if circulation_blade(i) <= 0;
%             circulation_blade(i) = 0;
%         end
%     end
    
    if timestep < 2
        TE_convected = TE_blade;
        TE_convected = TE_convected + repmat(dt*u_inf,1,length(TE_blade(1,:)),timestep);
        circulation_history = circulation_blade;
    else
        TE_convected = cat(3,TE_convected,TE_blade);
        TE_convected = TE_convected + repmat(dt*u_inf,1,length(TE_blade(1,:)),timestep);
        circulation_history = cat(3,circulation_history,circulation_blade);
    end
    
    B_ind = induction_factors(normalvectors_upd,cp_coords_blade,vort_end_blade,TE_blade,eps,1);
    [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,dT] = aero_coefficients(N,span,rho,chord,surface_area,chord,u_combined,u_inf,circulation_blade,B_ind,alpha,alphadata,Clspline,Cdspline,panel_width,cp_spanwise,omega);
    Thrust_convergence_ll(timestep) = sum(dT);
    
    
    time = dt*(timestep-1);
    rotation_angle = omega*time; % blade rotation angle w.r.t. positive y-axis [rad]
end

%% Calculating several forces & coefficients for the now-converged flow over the rotor blade
load('BEM_results_prandtl.mat');


B_ind = induction_factors(normalvectors_upd,cp_coords_blade,vort_end_blade,TE_blade,eps,1);
[dL,C_L,dD_i,C_D_i,alpha_inflow,alpha_i,alpha_eff,phi,C_N,C_T,w_ind,u_magnitude,dQ,dQ_cumulative,C_Q,C_Q_cumulative,dP,C_P,dT] = aero_coefficients(N,span,rho,chord,surface_area,chord,u_combined,u_inf,circulation_blade,B_ind,alpha,alphadata,Clspline,Cdspline,panel_width,cp_spanwise,omega);
P = omega*dQ_cumulative(end);
CP = P/(0.5*rho*(norm(u_inf)^3)*pi*span^2);
%% Remarks on convergence:
% 1) The circulation curve becomes smoother when less elements are used. I
% have no idea why this is, maybe the numerical round-off errors induced by
% matlab are at fault. When the number of elements is increased strange
% oscillations start occurring that do not cause divergence but are very
% persistent. Weird stuff man

%% Plotting module
figure(1)
plot(cp_mu,circulation_history(:,1,1),cp_mu,circulation_history(:,1,3),cp_mu,circulation_history(:,1,5),cp_mu,circulation_history(:,1,7),cp_mu,circulation_history(:,1,10),cp_mu,circulation_history(:,1,14),cp_mu,circulation_history(:,1,18),cp_mu,circulation_history(:,1,22),cp_mu,circulation_history(:,1,26),cp_mu,circulation_history(:,1,30),cp_mu,circulation_history(:,1,35),cp_mu,circulation_history(:,1,40),cp_mu,circulation_history(:,1,45),cp_mu,circulation_history(:,1,50),cp_mu,circulation_history(:,1,60),cp_mu,circulation_history(:,1,70),cp_mu,circulation_history(:,1,80),cp_mu,circulation_history(:,1,90));%,cp_mu,circulation_history(:,1,100))
grid minor
title('Convergence of circulation distribution')
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Circulation [m2/s]')

figure(2)
plot(cp_mu,circulation_history(:,1,timestep))
grid minor
title('Circulation distribution')
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Circulation [m2/s]')

figure(3)
plot(cp_mu,C_N,mu_local,Cx_opt)
grid minor
title('Normal force coefficient distribution')
legend('Lifting line','BEM','Location','South')
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Normal force coefficient C_N [-]')

figure(4)
plot(cp_mu,C_T,mu_local,C_t_opt)
grid minor
title('Thrust force coefficient distribution')
legend('Lifting line','BEM')
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Thrust force coefficient C_T [-]')

% figure(5)
% plot(cp_mu,rad2deg(alpha_eff))
% grid minor
% title('Effective angle of attack distribution')
% xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
% ylabel('Effective angle of attack [deg]')

figure(6)
plot(cp_mu,C_L)
grid minor
title('Lift coefficient distribution')
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Lift coefficient C_L [-]')

figure(7)
plot(cp_mu,C_D_i)
grid minor
title('Induced drag distribution')
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Induced drag coefficient C_D_i [-]')

figure(8)
plot(cp_mu,w_ind)
grid minor
title('Induced blade-normal velocity distribution')
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Induced blade-normal velocity [m/s]')

figure(9)
subplot(2,1,1)
plot(cp_mu,dQ_cumulative,mu_local,Torque_opt)
grid minor
title('Torque distribution')
legend('Lifting line','BEM','Location','North')
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Torque [N*m]')
subplot(2,1,2)
plot(cp_mu,C_Q_cumulative,mu_local,C_q_opt)
grid minor
title('Torque coefficient distribution')
legend('Lifting line','BEM','Location','North')
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Torque coefficient [-]')

figure(10)
subplot(2,1,1)
plot(cp_mu,rad2deg(phi),mu_local,phi_opt)
grid minor
title('Inflow angle')
legend('Lifting line','BEM')
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Inflow angle \phi [deg]')
subplot(2,1,2)
plot(cp_mu,rad2deg(alpha_eff),mu_local,AoA_opt)
grid minor
title('Angle of attack')
legend('Lifting line','BEM')
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Angle of attack \alpha [deg]')

figure(11)
plot(cp_mu,C_P,mu_local,Cp_opt)
grid minor
title('Power coefficient distribution')
legend('Lifting line','BEM','Location','North')
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Power coefficient [-]')

figure(12)
plot(linspace(1,max_timestep,max_timestep),Thrust_convergence_ll)
grid minor
title('Thrust convergence history')
xlabel('Iteration number [-]')
ylabel('Thrust [N]')

%% Special plots

figure(13)
plot(cp_mu_40,circulation_history_40(:,1,timestep),cp_mu_80,circulation_history_80(:,1,timestep),cp_mu_120,circulation_history_120(:,1,timestep),cp_mu_160,circulation_history_160(:,1,timestep))
grid minor
title('Circulation distribution for different numbers of elements')
legend('40 elements','80 elements','120 elements','Location','South')
xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
ylabel('Circulation [m2/s]')

figure(14)
plot(linspace(1,max_timestep,max_timestep),Thrust_convergence_40,linspace(1,max_timestep,max_timestep),Thrust_convergence_80,linspace(1,max_timestep,max_timestep),Thrust_convergence_120)
grid minor
title('Thrust convergence history')
legend('40 elements','80 elements','120 elements','Location','South')
xlabel('Iteration number [-]')
ylabel('Thrust [N]')