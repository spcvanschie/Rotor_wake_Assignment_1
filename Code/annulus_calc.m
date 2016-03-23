function [W,phi_all,AoA_all,Cx,Cy,a_new,a_tan_new]=annulus_calc(N,U_inf,r,omega,chordangle,Clspline,Cdspline,blade_solidity,B,mu_local,lambda,mu_min)
% Initial guesses for a and a_tangential
a = 0.01*ones(1,N);
a_tan = 0.01*ones(1,N);

a_underrelax = 0.1; % induction factor underrelaxation factor
a_tan_underrelax = 0.1; % tangential induction factor underrelaxation factor

a_new = zeros(1,N);
a_tan_new = zeros(1,N);
phi_all = zeros(1,N);
AoA_all = zeros(1,N);
for i = (1:N)
    a_1 = a(i);
    a_tan_1 = a_tan(i);
    
    run = 1; % iteration variable
    while run>0
        mu_current = mu_local(i) % console output to show current spanwise coordinate
        W = sqrt((U_inf.*(1-a_1)).^2+(r(i).*omega.*(1+a_tan_1)).^2);
        sinphi = (U_inf*(1-a_1))./W;
        phi = asind(sinphi);
        cosphi = (r(i).*omega.*(1+a_tan_1))./W;

        AoA = phi - chordangle(i); % angle of attack [deg]
        Cx = ppval(Clspline,AoA).*cosphi+ppval(Cdspline,AoA).*sinphi;
        Cy = ppval(Clspline,AoA).*sinphi-ppval(Cdspline,AoA).*cosphi;
        thrust_a_prelim = blade_solidity(i).*Cx./(4.*(sinphi.^2));
        torque_a_tan_prelim = blade_solidity(i).*Cy./(4.*sinphi.*cosphi);
        a_calc_prelim = thrust_a_prelim./(1+thrust_a_prelim);
        f_tip_1 = (2/pi)*acos(exp(-((B/2)*((1-mu_local(i))/mu_local(i))*sqrt(1+((lambda*mu_local(i))^2)/((1-a_calc_prelim)^2)))));
        f_root_1 = (2/pi)*acos(exp(-((B/2)*((mu_local(i)-mu_min)/mu_local(i))*sqrt(1+((lambda*mu_local(i))^2)/((1-a_calc_prelim)^2)))));
        f_1 = f_tip_1*f_root_1;
        a_calc = thrust_a_prelim*((1-a_calc_prelim)^2)/(f_1*(1-a_calc_prelim*f_1));
        a_tan_calc_prelim = torque_a_tan_prelim./(1-torque_a_tan_prelim);
        a_tan_calc = torque_a_tan_prelim*(1-a_calc)*(1+a_tan_calc_prelim)/(f_1*(1-a_calc*f_1));
        if abs(a_calc-a_1)<0.01*abs(a_calc) 
            if abs(a_tan_calc-a_tan_1)<0.01*abs(a_tan_calc)
                run = 0;
            end
        end
        a_1_new = a_1 + a_underrelax*(a_calc - a_1);
        a_1 = a_1_new;
        a_tan_1_new = a_tan_1 + a_tan_underrelax*(a_tan_calc - a_tan_1);
        a_tan_1 = a_tan_1_new;
    end
    a_new(i) = a_calc;
    a_tan_new(i) = a_tan_calc;
    phi_all(i) = phi;
    AoA_all(i) = AoA;
end