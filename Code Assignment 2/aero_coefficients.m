function [dL,C_L,dD_i,C_D_i,alpha_i,C_N,C_T,w_ind,u_magnitude] = aero_coefficients(rho,chordlength,surface_area,u_inflow,circulation,B_ind,alpha,Clspline,Cdspline)
u_magnitude = zeros(length(circulation),1);
for i = 1:length(circulation)
    u_magnitude(i) = norm(u_inflow(:,i));
end
w_ind = B_ind*circulation;

dL = zeros(length(circulation),1);
C_L = zeros(length(circulation),1);
dD_i = zeros(length(circulation),1);
C_D_i = zeros(length(circulation),1);
alpha_i = zeros(length(circulation),1);
alpha_eff = zeros(length(circulation),1);
C_N = zeros(length(circulation),1);
C_T = zeros(length(circulation),1);
for i = 1:length(circulation)
    alpha_i(i) = -rad2deg((w_ind(i)/u_magnitude(i)));
    alpha_eff(i) = -alpha_i(i) - alpha(i);
    dL(i) = rho*circulation(i)*u_magnitude(i)*(surface_area(i)/chordlength(i));
    C_L(i) = ppval(Clspline,alpha_eff(i)); %circulation(i)/(0.5*chordlength(i)*u_magnitude(i));
    dD_i(i) = -rho*w_ind(i)*circulation(i)*(surface_area(i)/chordlength(i));
    C_D_i(i) = ppval(Cdspline,alpha_eff(i)); %-circulation(i)/(0.5*chordlength(i)*w_ind(i));
    C_N(i) = cos(deg2rad(alpha(i)))*C_L(i) + sin(deg2rad(alpha(i)))*C_D_i(i);
    C_T(i) = -sin(deg2rad(alpha(i)))*C_L(i) + cos(deg2rad(alpha(i)))*C_D_i(i);

end
end