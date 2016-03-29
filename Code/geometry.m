function [r,R,B,mu_min,mu_local,twist,chordlength,chordangle,omega,blade_solidity]= geometry(N,pitch,lambda,U_inf,mu_min,twist_par,chordlength_par,mintwist_par,tip_par)
% Input parameters: N (number of annuli), pitch (blade pitch angle)

R = 50; % blade radius [m]
B = 3; % number of blades [-]

omega = lambda*U_inf/R; % Rotational speed [1/s]

mu_local = cosspace(mu_min+(1-mu_min)/N,1-(1-mu_min)/N,N); % r/R per annulus [-]
r = mu_local*R; % local annulus radius [m]

% adjustable design parameters for second part of assignment
% ----------------------------------------------------------
twist = (twist_par-mintwist_par)*(1-mu_local)+mintwist_par; % twist for each annulus [deg]
chordlength = (chordlength_par-tip_par)*(1-mu_local)+tip_par; % chord length for each annulus [m]
chordangle = twist+pitch; % chord angle (beta) for each annulus [deg]
blade_solidity = (B./(2*pi*mu_local)).*chordlength/R; % blade solidity (sigma_r) [-]
% ----------------------------------------------------------

end