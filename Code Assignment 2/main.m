clear variables;
%% Rotor & flow parameter declaration
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
n = 20; % number of collocation points [-]
dt = 10; % time step [s]
timestep = 0; % initial time step number [-]
rotation_angle = 0; % initial azimuthal blade rotation angle [rad], defined as 0 along the y-axis
eps = 10^-8; % vortex core epsilon

%% Calculation of initial collocation point coordinates
cp_init_y = cosspace(blade_root,span,n+2,0);
% Removal of root and tip points from CP y-locations
cp_init_y = cp_init_y(cp_init_y~=blade_root);
cp_init_y = cp_init_y(cp_init_y~=span);

cp_init_z = zeros(1,n);

%% Generation of initial vortex element endpoint coordinates
vort_end_y = endpoints(blade_root,span,cp_init_y);
vort_end_x = zeros(1,length(vort_end_y)); % the vortex elements are located in the yz-plane
vort_end_z = zeros(1,length(vort_end_y));
vort_end_init = [vort_end_x; vort_end_y; vort_end_z];
vort_end_spanwise = vort_end_y; % radial position of every vortex element endpoint

%% Calculating (constant) geometrical parameters for every discretised element
chord = chordlength(span,rootchord,tipchord,cp_init_y);
alpha = blade_angle(span,roottwist,tiptwist,bladepitch,cp_init_y);
surface_area = surfarea(chord,vort_end_y);
cp_spanwise = cp_init_y; % radial position of every collocation point

%% Adjustment of collocation point x-coordinates and trailing edge vortex corners
cp_init_x = 0.5*chord;
cp_init = [cp_init_x; cp_init_y; cp_init_z];

TE_x = 0.75*chordlength(span,rootchord,tipchord,vort_end_y);
TE = [TE_x; vort_end_y; vort_end_z];

%% Calculation of initial collocation point normal vectors
normalvectors = normalvect(alpha,rotation_angle);

u_inf_repmat = repmat(u_inf,1,length(normalvectors(1,:)));

%% Calculating initial influence coefficient matrix & freestream and blade rotation influence vector
A_ind_init = induction_factors(normalvectors,cp_init,vort_end_init,TE,eps);
u_rotor = rotationvect(omega,rotation_angle,cp_spanwise,normalvectors);
Q_inf = freestreamvect(normalvectors,u_inf_repmat,u_rotor);
%Q_rot = rotationvect(omega,rotation_angle,cp_spanwise,normalvectors);

%% Solving the initial linear system of equations for the vortex circulations
circulations = -A_ind_init\(Q_inf);

%% Time marching section for t > 0

max_timestep = 50;

Q_ind_start = zeros(length(chord),max_timestep);
for timestep = 0:max_timestep
    cp_coords_blade = bladerotation(rotation_angle,cp_spanwise,cp_init);
    vort_end_blade = bladerotation(rotation_angle,vort_end_spanwise,vort_end_init);
    TE_blade = bladerotation(rotation_angle,vort_end_spanwise,TE);
    normalvectors_upd = normalvect(alpha,rotation_angle);
    u_rotor_upd = rotationvect(omega,rotation_angle,cp_spanwise,normalvectors_upd);
    
    
    A_ind_upd = induction_factors(normalvectors_upd,cp_coords_blade,vort_end_blade,TE_blade,eps);
    Q_inf = freestreamvect(normalvectors_upd,u_inf_repmat,u_rotor_upd);
    
    if timestep < 1
        Q_ind = zeros(length(cp_spanwise),1);
    else
        [Q_ind, Q_coefs] = wakeinduction(Q_ind_start,cp_coords_blade,normalvectors_upd,TE_convected,TE_blade,circulation_history,eps);
    end
    
    circulation_blade = -A_ind_upd\(Q_inf + Q_ind);
    
    if timestep < 1
        TE_convected = TE_blade;
        TE_convected = TE_convected + repmat(dt*u_inf,1,length(TE_blade(1,:)),timestep + 1);
        circulation_history = circulation_blade;
    else
        TE_convected = cat(3,TE_convected,TE_blade);
        TE_convected = TE_convected + repmat(dt*u_inf,1,length(TE_blade(1,:)),timestep + 1);
        circulation_history = cat(3,circulation_history,circulation_blade);
    end
    %% Timestep break
    
    time = dt*timestep;
    rotation_angle = omega*time; % blade rotation angle w.r.t. positive y-axis [rad]
end

%% Plotting module
figure(1)
grid on
plot(cp_spanwise,circulation_history(:,1,1),cp_spanwise,circulation_history(:,1,3),cp_spanwise,circulation_history(:,1,5),cp_spanwise,circulation_history(:,1,7),cp_spanwise,circulation_history(:,1,10),cp_spanwise,circulation_history(:,1,14),cp_spanwise,circulation_history(:,1,18),cp_spanwise,circulation_history(:,1,22),cp_spanwise,circulation_history(:,1,26),cp_spanwise,circulation_history(:,1,30),cp_spanwise,circulation_history(:,1,35),cp_spanwise,circulation_history(:,1,40),cp_spanwise,circulation_history(:,1,45),cp_spanwise,circulation_history(:,1,50))