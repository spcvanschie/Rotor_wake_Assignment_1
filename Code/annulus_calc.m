function [r,W,phi,AoA]=annulus_calc(a,a_tan,U_inf,mu_local,r,R,omega,chordangle)
r = mu_local*R; % local blade radius [m]

W = sqrt((U_inf.*(1-a)).^2+(r.*omega.*(1+a_tan)).^2);
sinphi = (U_inf*(1-a))./W;
phi = asind(sinphi);
cosphi = (r.*omega.*(1+a_tan))./W;

AoA = phi - chordangle; % angle of attack for each annulus [deg]
end