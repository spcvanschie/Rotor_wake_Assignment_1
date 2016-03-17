function [r,W]=annulus_calc(a,a_tan,U_inf,mu,R,omega)
r = mu*R; % local blade radius [m]

W = sqrt((U_inf*(1-a))^2+(r*omega*(1+a_tan))^2); 
end