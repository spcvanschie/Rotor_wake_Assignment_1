function [r,R,B,mu_min,mu_local,twist,chordlength,chordangle,omega]= geometry(N,pitch,lambda,U_inf)
% Input parameters: N (number of annuli), pitch (blade pitch angle)

R = 50; % blade radius [m]
B = 3; % number of blades [-]
mu_min = 0.2; % spanwise start of blade [-]

omega = lambda*U_inf/R; % Rotational speed [1/s]

mu_local = linspace(mu_min,1,N); % r/R per annulus [-]
r = mu_local*R; % local annulus radius [m]
% adjustable design parameters for second part of assignment
% ----------------------------------------------------------
twist = 14*(1-mu_local); % twist for each annulus [deg]
chordlength = 3*(1-mu_local)+1; % chord length for each annulus [m]
chordangle = twist+pitch; % chord angle (beta) for each annulus [deg]
% ----------------------------------------------------------

end