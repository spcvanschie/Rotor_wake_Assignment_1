function [W,phi,AoA,Cx,Cy,a_new,a_tan_new]=annulus_calc(N,a,a_tan,U_inf,r,omega,chordangle,Clspline,Cdspline,blade_solidity)

a_new=zeros(1,N);
a_tan_new=zeros(1,N);
for i = (1:N)
    a_1 = a(i);
    a_tan_1 = a_tan(i);
    
    run = 1; % iteration variable
    while run>0
        W = sqrt((U_inf.*(1-a_1)).^2+(r(i).*omega.*(1+a_tan_1)).^2);
        sinphi = (U_inf*(1-a_1))./W;
        phi = asind(sinphi);
        cosphi = (r(i).*omega.*(1+a_tan_1))./W;

        AoA = phi - chordangle(i); % angle of attack [deg]
        Cx = ppval(Clspline,AoA).*cosphi+ppval(Cdspline,AoA).*sinphi;
        Cy = ppval(Clspline,AoA).*sinphi-ppval(Cdspline,AoA).*cosphi;
        thrust_a = blade_solidity(i).*Cx./(4.*sinphi.^2);
        torque_a_tan = blade_solidity(i).*Cy./(4.*sinphi.*cosphi);
        a_calc = thrust_a./(1+thrust_a);
        a_tan_calc = torque_a_tan./(1-torque_a_tan);
        if abs(a_calc-a_1)<0.01.*a_calc 
            if abs(a_tan_calc-a_tan_1)<0.01.*a_tan
                run = 0;
            end
        end
        a_1 = a_calc;
        a_tan_1 = a_tan_calc;
    end
    a_new(i) = a_calc;
    a_tan_new(i) = a_tan_calc;
end