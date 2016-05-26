function [dL,C_L,dD_i,C_D_i,alpha_i,C_N,C_T] = aero_coefficients(rho,chordlength,surface_area,u_inflow,circulation,B_ind)
u_magnitude = zeros(length(circulation),1);
for i = 1:length(circulation)
    u_magnitude = norm(u_inflow(:,i));
end
w_ind = B_ind*circulation;

dL = rho*circulation.*u_magnitude.*(surface_area./chordlength)';
C_L = circulation./(0.5*chordlength'.*u_magnitude);
dD_i = -rho*w_ind.*circulation.*(surface_area./chordlength)';
C_D_i = -circulation./(0.5*chordlength'.*w_ind);
alpha_i = -rad2deg(atan(w_ind./u_magnitude));
C_N = cos(deg2rad(alpha_i)).*C_L + sin(deg2rad(alpha_i)).*C_D_i;
C_T = -sin(deg2rad(alpha_i)).*C_L + cos(deg2rad(alpha_i)).*C_D_i;
end