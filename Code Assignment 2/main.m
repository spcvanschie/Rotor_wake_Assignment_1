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

%% Discretisation parameter declaration
n = 9; % number of collocation points [-]
dt = 0.1; % timestep [s]
rotation_angle = 0; % initial azimuthal blade rotation angle
eps = 10^-4; % vortex core epsilon

%% Calculation of initial collocation point coordinates
cp_init_y = cosspace(blade_root,span,n+2,1);
% Removal of root and tip points from CP y-locations
cp_init_y = cp_init_y(cp_init_y~=blade_root);
cp_init_y = cp_init_y(cp_init_y~=span);

cp_init_z = zeros(1,n);

%% Generation of vortex element endpoint coordinates
vort_end_y = endpoints(blade_root,span,cp_init_y);
vort_end_x = zeros(1,length(vort_end_y)); % the vortex elements are located in the yz-plane
vort_end_z = zeros(1,length(vort_end_y));
vort_end_init = [vort_end_x; vort_end_y; vort_end_z];

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
normalvect = normalvect(alpha,rotation_angle);

%% Calculating influence coefficient matrix & freestream influence vector
A_ind = induction_factors(normalvect,cp_init,vort_end_init,TE,eps);
Q_inf = freestreamvect(normalvect,u_inf);