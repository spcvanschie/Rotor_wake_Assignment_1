function [dL,C_L,dD_i,C_D_i,alpha_inflow,alpha_i,alpha_eff,phi,C_N,C_T,w_ind,u_magnitude,dQ,dQ_cumulative,C_Q,C_Q_cumulative,dP,C_P,dT] = aero_coefficients(N,R,rho,chordlength,surface_area,chord,u_inflow,u_freestream,circulation,B_ind,alpha,alphadata,Clspline,Cdspline,panel_width,cp_spanwise,omega)
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
alpha_inflow = zeros(length(circulation),1);
phi = zeros(length(circulation),1);
C_N = zeros(length(circulation),1);
C_T = zeros(length(circulation),1);
dQ = zeros(length(circulation),1);
dQ_cumulative = zeros(length(circulation),1);
C_Q = zeros(length(circulation),1);
C_Q_cumulative = zeros(length(circulation),1);
dT = zeros(length(circulation),1);
dP = zeros(length(circulation),1);
C_P = zeros(length(circulation),1);


for i = 1:length(circulation)
    phi(i) = atan(u_freestream(1)/(cp_spanwise(i)*omega));
    alpha_inflow(i) = phi(i) - deg2rad(alpha(i));
    alpha_i(i) = (-w_ind(i)/u_magnitude(i));
    alpha_eff(i) = alpha_inflow(i) - alpha_i(i);
    dL(i) = rho*circulation(i)*u_magnitude(i)*panel_width(i);
%     if alpha_eff(i) >= max(alphadata);
%         C_L(i) = ppval(Clspline,max(alphadata));
%     end
%     if rad2deg(alpha_eff(i)) <= min(alphadata);
%         C_L(i) = ppval(Clspline,min(alphadata));
%     end
%     if rad2deg(alpha_eff(i)) > min(alphadata) && rad2deg(alpha_eff(i)) < max(alphadata)
    C_L(i) = ppval(Clspline,rad2deg(alpha_eff(i))); %circulation(i)/(0.5*chordlength(i)*u_magnitude(i));
%     end
    dD_i(i) = -rho*w_ind(i)*circulation(i)*panel_width(i);
    C_D_i(i) = ppval(Cdspline,rad2deg(alpha_eff(i))); %-circulation(i)/(0.5*chordlength(i)*w_ind(i));
    C_N(i) = cos((phi(i)))*C_L(i) + sin((phi(i)))*C_D_i(i);
    C_T(i) = sin((phi(i)))*C_L(i) - cos((phi(i)))*C_D_i(i);
    dQ(i) = (dL(i)*sin(phi(i))-dD_i(i)*cos(phi(i)))*cp_spanwise(i);%0.5*rho*(u_magnitude(i)^2)*N*chord(i)*cp_spanwise(i)*(C_L(i)*sin((phi(i)))-C_D_i(i)*cos((phi(i))))*panel_width(i);
    dQ_cumulative(i) = sum(dQ(1:i));
    C_Q(i) = dQ(i)/(0.5*rho*(norm(u_freestream)^2)*pi*(R^3));
    C_Q_cumulative(i) = dQ_cumulative(i)/(0.5*rho*(norm(u_freestream)^2)*pi*(R^3));
    dT(i) = (dL(i)*cos(phi(i))+dD_i(i)*sin(phi(i)));%0.5*rho*(u_magnitude(i)^2)*N*chord(i)*(C_L(i)*cos((phi(i)))+C_D_i(i)*sin((phi(i))))*panel_width(i);
    dP(i) = dQ_cumulative(i)*omega;
    C_P(i) = dP(i)/(0.5*rho*(norm(u_freestream)^3)*pi*(R^2));
end
end