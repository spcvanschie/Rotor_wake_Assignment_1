%% Rotor & flow parameter declaration
span = 50; % span of the rotor blade [m]
blade_root = 0.2*span; % location of blade root [m]
rootchord = 4; % root chord length (at mu = 0) [m]
tipchord = 2; % tip chord length [m]
bladepitch = -1.5; % blade pitch angle [deg]
roottwist = 1; % root twist angle [deg]
tiptwist = 0; % tip twist angle [deg]
u_inf = 10; % freestream velocity [m/s]
tsr = 8; % Tip Speed Ratio [-]

%% Discretisation parameter declaration
n = 9; % number of collocation points [-]
dt = 0.1; % timestep [s]

%% Calculation of initial collocation point coordinates
cp_init_y = cosspace(blade_root,span,n+2,1);
% Removal of root and tip points from CP y-locations
cp_init_y = cp_init_y(cp_init_y~=blade_root);
cp_init_y = cp_init_y(cp_init_y~=span);

cp_init_x = zeros(1,n);
cp_init_z = zeros(1,n);
cp_init = [cp_init_x; cp_init_y; cp_init_z];

%% Generation of vortex element endpoint coordinates
vort_end_y = endpoints(blade_root,span,cp_init_y);
vort_end_x = zeros(1,length(vort_end_y));
vort_end_z = zeros(1,length(vort_end_y));
vort_end_init = [vort_end_x; vort_end_y; vort_end_z];

%% Calculating (constant) geometrical parameters for every discretised element
chord = chordlength(span,rootchord,tipchord,cp_init_y);
alpha = blade_angle(span,roottwist,tiptwist,bladepitch,cp_init_y);
surface_area = area(chord,vort_end_y)