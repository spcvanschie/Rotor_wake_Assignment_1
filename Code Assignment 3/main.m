clear variables;
%close all;

%% Declaration of variables
eps = 10^-6; % vortex core epsilon [m]
c = 1; % chord length [m]
alpha_amplitude = 5; % amplitude of angle of attack changes [deg]
alpha_0 = 3; % steady angle of attack [deg]
V_inf = [10; 0]; % freestream velocity components [m/s]

n = 5; % number of panels used [-]
dt = 0.0001; % time discretisation step [s]
max_timestep = 1000; % maximum number of time steps taken [-]
k = 0; % reduced frequency [-]

omega = k*(2*norm(V_inf))/c; % rotational frequency [rad/s]

alpha = zeros(1,max_timestep+1);
alphadot = zeros(1,max_timestep+1);

gamma_previous = zeros(1,n+1);
wake_coords = zeros(2,max_timestep+1);
vortex_strength = zeros(1,max_timestep+1);
CL = zeros(1,max_timestep+1);

for i= 0:max_timestep
    %% Computation of vortex and collocation point coordinates
    alpha(i+1) = alpha_amplitude*sin(i*dt*omega) + alpha_0;
    alphadot(i+1) = omega*alpha_amplitude*cos(i*dt*omega);
    [panel_TE,vort_chord,cp_chord,vort_coords,cp_coords,TE_coords] = coordinates(n,c,alpha(i+1));
    
    %% Computation of flat plate normal vector
    normal = normalvect(alpha(i+1));

    %% Computation of panel influence coefficients and circulations
    [a,RHS,gamma] = circ(i,n,V_inf,cp_coords,vort_coords,normal,alpha(i+1),alphadot(i+1),gamma_previous,wake_coords,eps);
    wake_coords(:,i+1) = TE_coords+0.25*dt*norm(V_inf)*[cos(deg2rad(alpha(i+1)));-sin(deg2rad(alpha(i+1)))];
    %% Update several quantities for the calculations during the next time step
    gamma_previous = gamma;
    vortex_strength(i+1) = gamma(end);
    CL(i+1) = sum(gamma(1:end-1))/(0.5*c*norm(V_inf));
    
    [wake_coords,wake_movement,wake_change,wake_convection] = wake_convect(n,i,dt,wake_coords,gamma,vortex_strength,vort_coords,V_inf,eps);
    
    %convect_matrix = cat(1,V_inf(1)*dt*ones(1,i+1),zeros(1,i+1));
    %wake_coords = wake_coords + cat(2,convect_matrix,zeros(2,max_timestep-i));
        
end

%% Calculating parameters at various points in the flow field
x_table = linspace(-c,2*c,50);
z_table = linspace(-c,2*c,50);
x_vel = zeros(length(z_table),length(x_table));
z_vel = zeros(length(z_table),length(x_table));
[x_vel,z_vel,pos_x,pos_z,vel_mag] = velocityfield(x_table,z_table,x_vel,z_vel,gamma,vortex_strength,vort_coords,wake_coords,eps);


%% Plots
figure(1)
plot(linspace(0,max_timestep,max_timestep+1),vortex_strength)
title('Strength of shed vortices')
xlabel('Time step [-]')

figure(2)
plot(cp_coords(1,:),cp_coords(2,:),'r+',wake_coords(1,:),wake_coords(2,:),'bo')
grid on

figure(3)
plot(cp_chord,gamma(1:end-1))

figure(4)
surf(pos_x,pos_z,vel_mag,'edgecolor','none')
view(2)
colormap(jet)
%alpha(0.1)
colorbar
caxis([0 10])

figure(5)
plot(linspace(0,max_timestep,max_timestep+1),CL)
title('Airfoil lift coefficient')
xlabel('Time step [-]')
ylabel('Airfoil lift coefficient C_l [1/m]')