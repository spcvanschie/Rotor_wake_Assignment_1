function [W,phi_all,AoA_all,Cx_all,Cy,a_new,a_tan_new,Q_all,Cq_all,Thrust_all,Cp_all,P,Thrust_convergence,Prandtl_all]=annulus_calc(rho,N,U_inf,r,R,omega,chordlength,chordangle,Clspline,Cdspline,blade_solidity,B,mu_local,lambda,mu_min,optimise,a_defined)
delta_mu = (1-mu_min-(2*(1-mu_min)/N))/N; % width of each annulus [-]

% Initial guesses for a and a_tangential
if optimise < 1
    a = 0.01*ones(1,N);
    a_underrelax = 0.1; % induction factor underrelaxation factor
else
    a = a_defined;
    a_underrelax = 0; % induction factor underrelaxation factor
end

a_tan = 0.01*ones(1,N);
a_tan_underrelax = 0.1; % tangential induction factor underrelaxation factor

a_new = zeros(1,N);
a_tan_new = zeros(1,N);
phi_all = zeros(1,N);
AoA_all = zeros(1,N);
Cx_all = zeros(1,N);
Thrust_all = zeros(1,N);
Q_all = zeros(1,N);
Cq_all = zeros(1,N);
Thrust_convergence = zeros(10,N);
Prandtl_all = zeros(1,N);
for i = (1:N)
    a_1 = a(i);
    a_tan_1 = a_tan(i);
    
    run = 1; % iteration variable
    iteration = 1; % iteration number
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
        if optimise < 1
            a_calc_prelim = thrust_a_prelim./(1+thrust_a_prelim);
        else
            a_calc_prelim = a_defined(i);
        end
        f_tip_1 = (2/pi)*acos(exp(-((B/2)*((1-mu_local(i))/mu_local(i))*sqrt(1+((lambda*mu_local(i))^2)/((1-a_calc_prelim)^2)))));
        f_root_1 = (2/pi)*acos(exp(-((B/2)*((mu_local(i)-mu_min)/mu_local(i))*sqrt(1+((lambda*mu_local(i))^2)/((1-a_calc_prelim)^2)))));
        f_1 = f_tip_1*f_root_1;
        if optimise < 1
            a_calc = thrust_a_prelim*((1-a_calc_prelim)^2)/(f_1*(1-a_calc_prelim*f_1));
            a_tan_calc_prelim = torque_a_tan_prelim./(1-torque_a_tan_prelim);
            a_tan_calc = torque_a_tan_prelim*(1-a_calc)*(1+a_tan_calc_prelim)/(f_1*(1-a_calc*f_1));
        else
            a_calc = a_calc_prelim;
            a_tan_calc_prelim = torque_a_tan_prelim./(1-torque_a_tan_prelim);
            a_tan_calc = torque_a_tan_prelim*(1-a_defined(i))*(1+a_tan_calc_prelim)/(f_1*(1-a_defined(i)*f_1));
        end
        if abs(a_calc-a_1)<0.01*abs(a_calc) 
            if abs(a_tan_calc-a_tan_1)<0.01*abs(a_tan_calc)
                run = 0;
            end
        end
        a_1_new = a_1 + a_underrelax*(a_calc - a_1);
        a_1 = a_1_new;
        a_tan_1_new = a_tan_1 + a_tan_underrelax*(a_tan_calc - a_tan_1);
        a_tan_1 = a_tan_1_new;
        Thrust_convergence(iteration,i) = 0.5*rho*(W^2)*B*chordlength(i)*Cx*R*delta_mu;
        iteration = iteration + 1;
    end
    % After the induction factors have been determined we can calculate the forces on each annulus
    
    Thrust_dmu = 0.5*rho*(W^2)*B*chordlength(i)*Cx*R*delta_mu;
    Torque_dmu = 0.5*rho*(W^2)*B*chordlength(i)*r(i)*Cy*R*delta_mu;
    Thrust_all(i) = Thrust_dmu;
    Q_all(i) = Torque_dmu;
    Cq_all(i) = Torque_dmu/(0.5*rho*(U_inf^2)*pi*(R^3));
    a_new(i) = a_calc;
    a_tan_new(i) = a_tan_calc;
    phi_all(i) = phi;
    AoA_all(i) = AoA;
    Cx_all(i) = Cx;
    Prandtl_all(i) = f_1;
end
Q_all_real = cumsum(Q_all);
Cq_all_real = cumsum(Cq_all);
Q_all = Q_all_real;
Cq_all = Cq_all_real;
Cp_all = Cq_all.*lambda;
P = Q_all(end)*omega;
end