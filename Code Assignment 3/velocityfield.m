function [x_vel,x_vel_ind,z_vel,z_vel_ind,pos_x,pos_z,vel_mag,CP] = velocityfield(x_table,z_table,x_vel,x_vel_ind,z_vel,z_vel_ind,gamma,vortex_strength,vort_coords,wake_coords,V_inf,eps)
pos_x = zeros(length(z_table),length(x_table));
pos_z = zeros(length(z_table),length(x_table));
vel_mag = zeros(length(z_table),length(x_table));
CP = zeros(length(z_table),length(x_table));
for i = 1:length(x_table)
    for k = 1:length(z_table)
        [x_vel_ind(k,i),x_vel(k,i),z_vel_ind(k,i),z_vel(k,i)] = velocity(x_table(i),z_table(k),gamma,vortex_strength,vort_coords,wake_coords,V_inf,eps);        
        pos_x(k,i) = x_table(i);
        pos_z(k,i) = z_table(k);
        vel_mag(k,i) = sqrt((x_vel(k,i))^2 + (z_vel(k,i))^2);
        CP(k,i) = 1 - (vel_mag(k,i)/norm(V_inf))^2;
    end
end
end