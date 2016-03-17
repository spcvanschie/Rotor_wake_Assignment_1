function [W,phi,AoA,Cx,Cy,a_new,a_tan_new]=annulus_calc(a,a_tan,U_inf,r,omega,chordangle,Clspline,Cdspline,blade_solidity)
run = 1; % iteration variable

for i
a_1 = a;
a_tan_1 = a_tan;

while run>0
    W = sqrt((U_inf.*(1-a_1)).^2+(r.*omega.*(1+a_tan_1)).^2);
    sinphi = (U_inf*(1-a_1))./W;
    phi = asind(sinphi);
    cosphi = (r.*omega.*(1+a_tan_1))./W;

    AoA = phi - chordangle; % angle of attack for each annulus [deg]
    Cx = ppval(Clspline,AoA).*cosphi+ppval(Cdspline,AoA).*sinphi;
    Cy = ppval(Clspline,AoA).*sinphi-ppval(Cdspline,AoA).*cosphi;
    thrust_a = blade_solidity.*Cx./(4.*sinphi.^2);
    torque_a_tan = blade_solidity.*Cy./(4.*sinphi.*cosphi);
    a_new = thrust_a./(1+thrust_a);
    a_tan_new = torque_a_tan./(1-torque_a_tan);
    if abs(a_new-a_1)<0.01.*a_new && abs(a_tan_new-a_tan_1)<0.01.*a_tan
        run = 0;
    end
    a_1 = a_new;
    a_tan_1 = a_tan_new;
end
end