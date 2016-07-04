clear variables;
%close all;

%% Declaration of variables
eps = 10^-6; % vortex core epsilon [m]
c = 1; % chord length [m]
alpha_amplitude = 7; % amplitude of angle of attack changes [deg]
alpha_0 = 7; % steady angle of attack [deg]
V_inf = [10; 0]; % freestream velocity components [m/s]

n = 30; % number of panels used [-]
dt = 0.05; % time discretisation step [s]
max_timestep = 500; % maximum number of time steps taken [-]
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
    wake_coords(:,i+1) = TE_coords+0.1*dt*norm(V_inf)*[cos(deg2rad(alpha(i+1)));-sin(deg2rad(alpha(i+1)))];
    %% Update several quantities for the calculations during the next time step
    gamma_previous = gamma;
    vortex_strength(i+1) = gamma(end);
    CL(i+1) = sum(gamma(1:end-1))/(0.5*c*norm(V_inf));
    
    [wake_coords,wake_movement,wake_change,wake_convection] = wake_convect(n,i,dt,wake_coords,gamma,vortex_strength,vort_coords,V_inf,eps);
    
    %convect_matrix = cat(1,V_inf(1)*dt*ones(1,i+1),zeros(1,i+1));
    %wake_coords = wake_coords + cat(2,convect_matrix,zeros(2,max_timestep-i));
        
end

%% Calculating parameters at various points in the flow field
x_table = linspace(-c,2*c,100);
z_table = linspace(-1.5*c,1.5*c,100);
x_vel = zeros(length(z_table),length(x_table));
x_vel_ind = zeros(length(z_table),length(x_table));
z_vel = zeros(length(z_table),length(x_table));
z_vel_ind = zeros(length(z_table),length(x_table));
[x_vel,x_vel_ind,z_vel,z_vel_ind,pos_x,pos_z,vel_mag,CP] = velocityfield(x_table,z_table,x_vel,x_vel_ind,z_vel,z_vel_ind,gamma,vortex_strength,vort_coords,wake_coords,V_inf,eps);
[x_mesh_pos, z_mesh_pos] = meshgrid(x_table,z_table);

%% Plots
figure(1)
plot(linspace(0,max_timestep,max_timestep+1),vortex_strength)
title('Strength of shed vortices')
xlabel('Time step [-]')

figure(2)
plot(cp_coords(1,:),cp_coords(2,:),'r+',wake_coords(1,:),wake_coords(2,:),'bo')
title('Wake vortex pattern')
legend('Control points','Shed vortex cores')
xlabel('x [m]')
ylabel('y [m]')
grid on

figure(3)
plot(cp_chord,gamma(1:end-1))
title('Circulation distribution')
xlabel('\frac{x}{c} [-]')
ylabel('Circulation [m^2/s]')

figure(4)
surf(pos_x,pos_z,vel_mag,'edgecolor','none')
hold on
plot3(cp_coords(1,:),cp_coords(2,:),20*ones(length(cp_coords(1,:))),'r+')
title('Velocity magnitude distribution around flat plate')
view(2)
colormap(jet)
%alpha(0.1)
bar = colorbar;
bar.Label.String = 'Velocity magnitude [m/s]';
caxis([8 12])
xlabel('x [m]')
ylabel('y [m]')
hold off

figure(5)
quiver(x_mesh_pos,z_mesh_pos,x_vel_ind,z_vel_ind)
hold on
%plot(cp_coords(1,:),cp_coords(2,:),20*ones(length(cp_coords(1,:))),'r+')
%view(2)
%colormap(jet)
%alpha(0.1)
%colorbar
%caxis([9 11])
xlabel('x [m]')
ylabel('y [m]')
hold off

figure(6)
plot(linspace(0,max_timestep,max_timestep+1),CL)
title('Flat plate lift coefficient')
xlabel('Time step [-]')
ylabel('Flat plate lift coefficient C_l [1/m]')

figure(7)
surf(pos_x,pos_z,CP,'edgecolor','none')
hold on
plot3(cp_coords(1,:),cp_coords(2,:),ones(length(cp_coords(1,:))),'r+')
title('Pressure coefficient distribution around flat plate')
view(2)
colormap(jet)
%alpha(0.1)
bar2 = colorbar;
bar2.Label.String = 'Pressure coefficient C_P [-]';
caxis([-0.5 0.5])
hold off