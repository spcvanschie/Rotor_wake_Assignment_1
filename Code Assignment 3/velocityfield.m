function [x_vel,z_vel,pos_x,pos_z,vel_mag] = velocityfield(x_table,z_table,x_vel,z_vel,gamma,vortex_strength,vort_coords,wake_coords,eps)
pos_x = zeros(length(z_table),length(x_table));
pos_z = zeros(length(z_table),length(x_table));
vel_mag = zeros(length(z_table),length(x_table));
for i = 1:length(x_table)
    for k = 1:length(z_table)
        [x_vel(k,i),z_vel(k,i)] = velocity(x_table(i),z_table(k),gamma,vortex_strength,vort_coords,wake_coords,eps);        
        pos_x(k,i) = x_table(i);
        pos_z(k,i) = z_table(k);
        vel_mag(k,i) = sqrt((x_vel(k,i))^2 + (z_vel(k,i))^2);
    end
end
end