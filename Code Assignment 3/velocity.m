function [x_vel_ind,z_vel_ind] = velocity(x,z,plate_circulation,wake_circulation,vort_coords,wake_coords,V_inf,eps)
v_ind = [0;0];
for i = 1:length(vort_coords(1,:))
    v_ind = v_ind + VOR2D(plate_circulation(i),x,z,vort_coords(1,i),vort_coords(2,i),eps);    
end
for j = 1:length(wake_circulation)
    v_ind = v_ind + VOR2D(wake_circulation(j),x,z,wake_coords(1,j),wake_coords(2,j),eps);
end
x_vel_ind = v_ind(1) + V_inf(1);
z_vel_ind = v_ind(2) + V_inf(2);
end