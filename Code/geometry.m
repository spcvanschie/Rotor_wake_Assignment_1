function [R,B,mu_min,mu,twist,chordlength]= geometry(N,pitch)
% Input parameters: N (number of annuli), pitch (blade pitch angle)

R = 50; % blade radius [m]
B = 3; % number of blades [-]
mu_min = 0.2; % spanwise start of blade [-]

mu = linspace(mu_min,1,N); % r/R per annulus [-]
twist = 14*(1-mu); % twist per annulus [deg]
chordlength = 3*(1-mu)+1; % chord length per annulus [m]

end