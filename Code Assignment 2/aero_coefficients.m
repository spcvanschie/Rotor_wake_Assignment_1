function [dL,C_L,dD_i,C_D_i,alpha_i,C_N,C_T,w_ind,u_magnitude] = aero_coefficients(rho,chordlength,surface_area,u_inflow,circulation,B_ind)
u_magnitude = zeros(length(circulation),1);
for i = 1:length(circulation)
    u_magnitude(i) = norm(u_inflow(:,i));
end
w_ind = B_ind'*circulation;

dL = zeros(length(circulation),1);
C_L = zeros(length(circulation),1);
dD_i = zeros(length(circulation),1);
C_D_i = zeros(length(circulation),1);
alpha_i = zeros(length(circulation),1);

for i = 1:length(circulation)
    dL(i) = rho*circulation(i)*u_magnitude(i)*(surface_area(i)/chordlength(i));
    C_L(i) = circulation(i)/(0.5*chordlength(i)*u_magnitude(i));
    dD_i(i) = -rho*w_ind(i)*circulation(i)*(surface_area(i)/chordlength(i));
    C_D_i(i) = -circulation(i)/(0.5*chordlength(i)*w_ind(i));
    alpha_i(i) = -rad2deg((w_ind(i)/u_magnitude(i)));
end
C_N = cos(deg2rad(alpha_i)).*C_L + sin(deg2rad(alpha_i)).*C_D_i;
C_T = -sin(deg2rad(alpha_i)).*C_L + cos(deg2rad(alpha_i)).*C_D_i;
end