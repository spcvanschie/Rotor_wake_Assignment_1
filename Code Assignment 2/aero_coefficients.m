function [L,C_L,D_i,C_D_i] = aero_coefficients(rho,chordlength,u_inflow,circulation,B_ind)
for i = 1:length(circulation)
    u_magnitude = norm(u_inflow(:,i));
end
L = - rho*circulation.*u_magnitude
C_L = circulation./(0.5*chordlength.*u_magnitude);
D_i =
C_D_i = 

end