clear all;

N = [100]; % number of annuli
lambda = [6,8,10]; % tip speed ratio [-]
U_inf = 10; % freestream velocity [m/s]
mu_min = 0.2; % spanwise start of blade [-]
mu_min_hub = mu_min; % spanwise edge of rotor hub [-]
rho = 1.225; % air density [kg/m^3]

% Optimisation parameters
pitch = 2; % blade pitch angle [deg]
twist_par = 14; % twist for each annulus [deg]
mintwist_par = 0; % tip twist
chordlength_par = 3; % chord length for each annulus [m]
tip_par = 1; % tip chord length

baseline = 0; % bogey stagement to either use or ignore the baseline studycase part of the code
optimise = 1; % bogey statement to either use or ignore the optimisation part of the code

if baseline > 0
    mu=cell(length(N),3);
    a_all=cell(length(N),3);
    a_tan_all=cell(length(N),3);
    AoA_all=cell(length(N),3);
    phi_all=cell(length(N),3);
    C_t_all=cell(length(N),3);
    C_n_all=cell(length(N),3);
    Thrust_all = cell(length(N),3);
    Torque_all = cell(length(N),3);
    Cq_all = cell(length(N),3);
    Power_baseline = ones(1,3);
    Cp_baseline = ones(1,3);
    Thrust_convergence = cell(length(N),3);
    Prandtl_all = cell(length(N),3);
    for j = (1:length(N))
        for i = (1:3)
            % calculate interpolation splines for airfoil Cl and Cd
            [Clspline,Cdspline,alphadata,Cldata,Cddata]=airfoil_liftdrag();

            % calculate geometrical parameters for each annulus
            [r,R,B,mu_min,mu_local,twist,chordlength,chordangle,omega,blade_solidity_baseline]= geometry(N(j),pitch,lambda(i),U_inf,mu_min,twist_par,chordlength_par,mintwist_par,tip_par);

            % calculate annulus characteristics, contains iteration loop for induction factors
            [W,phi,AoA,Cx,Cy,a_new,a_tan_new,Torque,C_torque,Thrust,Cp_all,P,Thrust_converge,Prandtl]=annulus_calc(rho,N(j),U_inf,r,R,omega,chordlength,chordangle,Clspline,Cdspline,blade_solidity_baseline,B,mu_local,lambda(i),mu_min,mu_min_hub,0,0);

        % writing all obtained variables to the memory    
        mu{j,i} = mu_local;
        a_all{j,i} = a_new;
        a_tan_all{j,i} = a_tan_new;
        AoA_all{j,i} = AoA;
        phi_all{j,i} = phi;
        C_t_all{j,i} = Glauert(a_new);
        C_n_all{j,i} = Cx;
        Thrust_all{j,i} = Thrust;
        Torque_all{j,i} = Torque;
        Cq_all{j,i} = C_torque;
        Prandtl_all{j,i} = Prandtl;
        Power_baseline(1,i) = P;
        Cp_baseline(1,i) = P/(0.5*rho*(U_inf^3)*pi*R*R);
        Thrust_convergence{j,i} = Thrust_convergence_history(Thrust_converge,Thrust);
        end
    end
end

%% Blade optimisation section of code
pitch = 2; % blade pitch angle [deg]
lambda_optimise = 8;

% starting values of optimisation parameters
maxtwist_min = -20; % root twist angle [deg]
maxtwist_max = 20; % root twist angle [deg]
maxtwist_samples = 10; % number of samples for maxtwist
mintwist_min = 0;
mintwist_max = 0;
mintwist_samples = 1;
rootminustip_min = 1; % root chord length minus tip chord [m]
rootminustip_max = 2; % root chord length minus tip chord [m]
rootminustip_samples = 2; % number of samples for rootminustip
tip_min = 1;
tip_max = 2;
tip_samples = 2;
pitch_min = -15; % min pitch angle [deg]
pitch_max = 15; % max pitch angle [deg]
pitch_samples = 15; % number of pitch angle samples


maxtwist_range = linspace(maxtwist_min,maxtwist_max,maxtwist_samples);
mintwist_range = linspace(mintwist_min,mintwist_max,mintwist_samples);
rootminustip_range = linspace(rootminustip_min,rootminustip_max,rootminustip_samples);
tip_range = linspace(tip_min,tip_max,tip_samples);
pitch_range = linspace(pitch_min,pitch_max,pitch_samples);

if optimise > 0
    C_t_factor = 1;
    
    Power_data = zeros(length(maxtwist_range),length(rootminustip_range),length(mintwist_range),length(tip_range));
    maxtwist_data = zeros(length(maxtwist_range),length(rootminustip_range),length(mintwist_range),length(tip_range));
    mintwist_data = zeros(length(maxtwist_range),length(rootminustip_range),length(mintwist_range),length(tip_range));
    rootminustip_data = zeros(length(maxtwist_range),length(rootminustip_range),length(mintwist_range),length(tip_range));
    tip_data = zeros(length(maxtwist_range),length(rootminustip_range),length(mintwist_range),length(tip_range));
    for i = (1:length(maxtwist_range))
        for j = (1:length(rootminustip_range))
            for k = (1:length(mintwist_range))
                for l = (1:length(tip_range))
                
                    twist_par = maxtwist_range(i); % twist for each annulus [deg]
                    chordlength_par = rootminustip_range(j); % chord length for each annulus [m]
                    mintwist_par = mintwist_range(k);
                    tip_par = tip_range(l);
                    
                    run = 1;
                    C_t_underrelax = 0.5;
                    C_t_design = 0.75;
                    while run>0;
                        a = a_from_C_t(C_t_design)*ones(1,max(N));  
                        % calculate interpolation splines for airfoil Cl and Cd
                        [Clspline,Cdspline,alphadata,Cldata,Cddata]=airfoil_liftdrag();

                        % calculate geometrical parameters for each annulus
                        [r,R,B,mu_min,mu_local,twist,chordlength,chordangle,omega_opt,blade_solidity]= geometry(max(N),pitch,lambda_optimise,U_inf,mu_min,twist_par,chordlength_par,mintwist_par,tip_par);

                        % calculate annulus characteristics, contains iteration loop for induction factors
                        [W,phi,AoA,Cx,Cy,a_new,a_tan_new,Torque,C_torque,Thrust,Cp_opt,Power_optloop,Thrust_convergence_design]=annulus_calc(rho,max(N),U_inf,r,R,omega_opt,chordlength,chordangle,Clspline,Cdspline,blade_solidity,B,mu_local,lambda_optimise,mu_min,mu_min_hub,optimise,a);
                        
                        C_t_actual = mean(Glauert(a_new));
                        C_t_design = C_t_design +C_t_underrelax*(0.75 - C_t_actual);
                        if abs(mean(Glauert(a_new))-0.75) < 0.01
                            run = 0;
                        end
                    end
                    Power_data(i,j,k,l) = Power_optloop;
                    maxtwist_data(i,j,k,l) = maxtwist_range(i);
                    rootminustip_data(i,j,k,l) = rootminustip_range(j);
                    mintwist_data(i,j,k,l) = mintwist_range(k);
                    tip_data(i,j,k,l) = tip_range(l);
                end
            end
        end
    end
    
    %Power_interp = griddata(maxtwist_data,rootminustip_data,Power_data,maxtwist_data,rootminustip_data,'linear');
    
    % Optimize pitch angle
    [maxpower, index] = max(Power_data(:));
    [opt_maxtwistindex,opt_rootminustipindex,opt_mintwistindex,opt_tipindex] = ind2sub(size(Power_data),index);
    opt_maxtwist = maxtwist_range(opt_maxtwistindex);
    opt_rootminustip = rootminustip_range(opt_rootminustipindex);
    opt_mintwist = mintwist_range(opt_mintwistindex);
    opt_tip = tip_range(opt_tipindex);
    
    twist_par_opt = opt_maxtwist; % twist for each annulus [deg]
    chordlength_par_opt = opt_rootminustip; % chord length for each annulus [m]
    mintwist_par_opt = opt_mintwist;
    tip_par_opt = opt_tip;
    
    Power_pitch = zeros(1,length(pitch_range));
    for i = (1:length(pitch_range))
        % calculate geometrical parameters for each annulus
        [r,R,B,mu_min,mu_local,twist,chordlength,chordangle,omega_pitch,blade_solidity]= geometry(max(N),pitch_range(i),lambda_optimise,U_inf,mu_min,twist_par_opt,chordlength_par_opt,mintwist_par_opt,tip_par_opt);

        % calculate annulus characteristics, contains iteration loop for induction factors
        [W,phi,AoA,Cx,Cy,a_new,a_tan_new,Torque,C_torque,Thrust,Cp_opt,P,Thrust_convergence_design]=annulus_calc(rho,max(N),U_inf,r,R,omega_pitch,chordlength,chordangle,Clspline,Cdspline,blade_solidity,B,mu_local,lambda_optimise,mu_min,mu_min_hub,optimise,a);
        Power_pitch(i) = P;
    end
    [maxpower_pitch,index] = max(Power_pitch);
    optdesign_pitch = ind2sub(size(Power_pitch),index);
    opt_pitch = pitch_range(optdesign_pitch);
    
    % calculate geometrical parameters for each annulus
    [r,R,B,mu_min,mu_local,twist,chordlength,chordangle,omega_pitch,blade_solidity]= geometry(max(N),opt_pitch,lambda_optimise,U_inf,mu_min,twist_par_opt,chordlength_par_opt,mintwist_par_opt,tip_par_opt);

    % calculate annulus characteristics, contains iteration loop for induction factors
    [W,phi_opt,AoA_opt,Cx,Cy,a_opt,a_tan_opt,Torque,C_q_opt,Thrust,Cp_opt,Maxpower,Thrust_convergence_design,Prandtl]=annulus_calc(rho,max(N),U_inf,r,R,omega_pitch,chordlength,chordangle,Clspline,Cdspline,blade_solidity,B,mu_local,lambda_optimise,mu_min,mu_min_hub,optimise,a);
    Cp_maxpower = Maxpower/(0.5*rho*(U_inf^3)*pi*(R^2))
    C_t_opt = Glauert(a_new);    
end


%% Plotting section of code
% figure axis ranges
axis_alpha = [0.2 1 -2 15];
axis_phi = [0.2 1 -2 30];
axis_a = [0.2 1 0 1];
axis_a_tan = [0.2 1 0 0.1];
axis_C_t = [0.2 1 0 1.2];
axis_C_n = [0.2 1 0 1.2];
axis_C_q = [0.2 1 0 1.5];
axis_Cl = [-16 30 -1 1.5];
axis_Cd = [-16 30 0 0.7];

if baseline > 0
    figure(1)
    subplot(2,1,1)
    plot(mu_local,AoA_all{length(N),1},mu_local,AoA_all{length(N),2},mu_local,AoA_all{length(N),3})
    grid on
    title('Angle of Attack')
    legend('\lambda = 6','\lambda = 8','\lambda = 10')
    axis(axis_alpha)
    xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
    ylabel('Angle of Attack \alpha [deg]')
    subplot(2,1,2)
    plot(mu_local,phi_all{length(N),1},mu_local,phi_all{length(N),2},mu_local,phi_all{length(N),3})
    grid on
    title('Inflow angle')
    legend('\lambda = 6','\lambda = 8','\lambda = 10')
    axis(axis_phi)
    xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
    ylabel('Inflow angle \phi [deg]')
    
    figure(2)
    plot(mu_local,Prandtl_all{length(N),1},mu_local,Prandtl_all{length(N),2},mu_local,Prandtl_all{length(N),3})
    grid on
    title('Combined tip and root corrections')
    legend('\lambda = 6','\lambda = 8','\lambda = 10')
    %axis(axis_C_t)
    xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
    ylabel('Correction factor [-]')
  
    figure(3)
    subplot(2,1,1)
    plot(mu_local,a_all{length(N),1},mu_local,a_all{length(N),2},mu_local,a_all{length(N),3})
    grid on
    title('Induction factor')
    legend('\lambda = 6','\lambda = 8','\lambda = 10')
    axis(axis_a)
    xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
    ylabel('Induction factor a [-]')
    subplot(2,1,2)
    plot(mu_local,a_tan_all{length(N),1},mu_local,a_tan_all{length(N),2},mu_local,a_tan_all{length(N),3})
    grid on
    title('Tangential induction factor')
    legend('\lambda = 6','\lambda = 8','\lambda = 10')
    axis(axis_a_tan)
    xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
    ylabel('Induction factor a'' [-]')

    figure(4)
    subplot(2,1,1)
    plot(mu_local,C_t_all{length(N),1},mu_local,C_t_all{length(N),2},mu_local,C_t_all{length(N),3})
    grid on
    title('Thrust coefficient')
    legend('\lambda = 6','\lambda = 8','\lambda = 10')
    axis(axis_C_t)
    xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
    ylabel('Thrust coefficient C_{T} [-]')
    subplot(2,1,2)
    plot(mu_local,C_n_all{length(N),1},mu_local,C_n_all{length(N),2},mu_local,C_n_all{length(N),3})
    grid on
    title('Normal force coefficient')
    legend('\lambda = 6','\lambda = 8','\lambda = 10')
    axis(axis_C_n)
    xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
    ylabel('Normal force coefficient C_{n} [-]')

    figure(5)
    subplot(2,1,1)
    plot(mu_local,Torque_all{length(N),1},mu_local,Torque_all{length(N),2},mu_local,Torque_all{length(N),3})
    grid on
    title('Torque')
    legend('\lambda = 6','\lambda = 8','\lambda = 10')
    %axis(axis_C_t)
    xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
    ylabel('Torque [N*m]')
    subplot(2,1,2)
    plot(mu_local,Cq_all{length(N),1},mu_local,Cq_all{length(N),2},mu_local,Cq_all{length(N),3})
    grid on
    title('Torque coefficient')
    legend('\lambda = 6','\lambda = 8','\lambda = 10')
    %axis(axis_C_q)
    xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
    ylabel('Torque coefficient C_{q} [-]')

    figure(6)
    plot(linspace(1,length(Thrust_convergence{length(N),1}),length(Thrust_convergence{length(N),1})),Thrust_convergence{length(N),1},linspace(1,length(Thrust_convergence{length(N),2}),length(Thrust_convergence{length(N),2})),Thrust_convergence{length(N),2},linspace(1,length(Thrust_convergence{length(N),3}),length(Thrust_convergence{length(N),3})),Thrust_convergence{length(N),3})
    grid on
    title('Thrust convergence history')
    legend('\lambda = 6','\lambda = 8','\lambda = 10')
    %axis(axis_C_t)
    xlabel('Iteration number')
    ylabel('Blade thrust [N]')

    figure(7)
    subplot(2,1,1)
    alpha_interpolation=linspace(min(alphadata),max(alphadata),500);
    plot(alphadata,Cldata,'o',alpha_interpolation,ppval(Clspline,alpha_interpolation))
    grid on
    title('C_{l} - \alpha data')
    legend('C_{l}-values','Interpolation data')
    axis(axis_Cl)
    xlabel('Angle of attack \alpha [deg]')
    ylabel('C_{l} [-]')
    subplot(2,1,2)
    plot(alphadata,Cddata,'o',alpha_interpolation,ppval(Cdspline,alpha_interpolation))
    grid on
    title('C_{d} - \alpha data')
    legend('C_{d}-values','Interpolation curve')
    axis(axis_Cd)
    xlabel('Angle of attack \alpha [deg]')
    ylabel('C_{d} [-]')
end

if optimise > 0
    if ndims(Power_data) < 3
        figure(8)
        surf(maxtwist_data,rootminustip_data,Power_data)
        xlabel('Maximum twist angle [deg]')
        ylabel('Chord length difference between root and tip [m]')
        zlabel('Power generated by rotor [W]')
    end
    
    figure(9)
    plot(pitch_range,Power_pitch)
    grid on
    title('Power generated for varying pitch angles')
    %axis(axis_Cd)
    xlabel('Pitch angle [deg]')
    ylabel('Power [W]')
    
    figure(10)
    subplot(2,1,1)
    plot(mu_local,AoA_opt)
    grid on
    title('Angle of attack of optimum rotor')
    xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
    ylabel('Angle of attack \alpha [deg]')
    subplot(2,1,2)
    plot(mu_local,phi_opt)
    grid on
    title('Inflow angle of optimum rotor')
    xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
    ylabel('Inflow angle \phi [deg]')
    
    figure(11)
    subplot(2,1,1)
    plot(mu_local,a_opt)
    grid on
    title('Induction factor of optimum rotor')
    xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
    ylabel('Induction factor [-]')
    subplot(2,1,2)
    plot(mu_local,a_tan_opt)
    grid on
    title('Axial induction factor of optimum rotor')
    xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
    ylabel('Axial induction factor [-]')
    
    figure(12)
    subplot(2,1,1)
    plot(mu_local,C_t_opt)
    grid on
    title('Thrust coefficient of optimum rotor')
    xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
    ylabel('Thrust coefficient [-]')
    subplot(2,1,2)
    plot(mu_local,C_q_opt)
    grid on
    title('Torque coefficient of optimum rotor')
    xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
    ylabel('Torque coefficient [-]')
end

if length(N)>1
    figure(13)
    plot(mu{1,2},a_all{1,2},mu{2,2},a_all{2,2},mu{3,2},a_all{3,2})
    grid on
    title('Induction factor (\lambda = 8)')
    legend('N = 100','N = 300','N = 500')
    axis(axis_a)
    xlabel('$\frac{r}{R} [-]$','Interpreter','LaTex')
    ylabel('Induction factor a [-]')
end