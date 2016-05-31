function [dL,C_L,dD_i,C_D_i,alpha_i,C_N,C_T,w_ind,u_magnitude,dQ,dQ_cumulative,C_Q,C_Q_cumulative,dT] = aero_coefficients(N,R,rho,chordlength,surface_area,chord,u_inflow,u_freestream,circulation,B_ind,alpha,Clspline,Cdspline,panel_width,cp_spanwise)
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
dQ = zeros(length(circulation),1);
dQ_cumulative = zeros(length(circulation),1);
C_Q = zeros(length(circulation),1);
C_Q_cumulative = zeros(length(circulation),1);
dT = zeros(length(circulation),1);
for i = 1:length(circulation)
    alpha_i(i) = -rad2deg((w_ind(i)/u_magnitude(i)));
    alpha_eff(i) = -alpha_i(i) - alpha(i);
    dL(i) = rho*circulation(i)*u_magnitude(i)*panel_width(i);
    C_L(i) = ppval(Clspline,alpha_eff(i)); %circulation(i)/(0.5*chordlength(i)*u_magnitude(i));
    dD_i(i) = -rho*w_ind(i)*circulation(i)*panel_width(i);
    C_D_i(i) = ppval(Cdspline,alpha_eff(i)); %-circulation(i)/(0.5*chordlength(i)*w_ind(i));
    C_N(i) = cos(deg2rad(alpha(i)))*C_L(i) + sin(deg2rad(alpha(i)))*C_D_i(i);
    C_T(i) = -sin(deg2rad(alpha(i)))*C_L(i) + cos(deg2rad(alpha(i)))*C_D_i(i);
    dQ(i) = 0.5*rho*(u_magnitude(i)^2)*N*chord(i)*cp_spanwise(i)*(dL(i)*sin(deg2rad(alpha_eff(i)))-dD_i(i)*cos(deg2rad(alpha_eff(i))))*panel_width(i);
    dQ_cumulative(i) = sum(dQ(1:i));
    C_Q(i) = dQ(i)/(0.5*rho*(norm(u_freestream)^2)*pi*(R^3));
    C_Q_cumulative(i) = dQ_cumulative(i)/(0.5*rho*(norm(u_freestream)^2)*pi*(R^3));
    dT(i) = 0.5*rho*(u_magnitude(i)^2)*N*chord(i)*(C_L(i)*cos(deg2rad(alpha_eff(i)))+C_D_i(i)*sin(deg2rad(alpha_eff(i))))*panel_width(i);
end
end