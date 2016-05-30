function [a_inf] = HSHOE(vect_norm,x_cp,x_a,x_b,x_c,x_d,circ,eps,drag);
[u_x1,u_y1,u_z1] = VORTXL(x_cp,x_a,x_b,circ,eps);
[u_x2,u_y2,u_z2] = VORTXL(x_cp,x_b,x_c,circ,eps);
[u_x3,u_y3,u_z3] = VORTXL(x_cp,x_c,x_d,circ,eps);
[u_x4,u_y4,u_z4] = VORTXL(x_cp,x_d,x_a,circ,eps);
if drag > 0
    u_x2 = 0;
    u_y2 = 0;
    u_z2 = 0;
    u_x4 = 0;
    u_y4 = 0;
    u_z4 = 0;
end
u_ind_x = u_x1 + u_x2 + u_x3 + u_x4;
u_ind_y = u_y1 + u_y2 + u_y3 + u_y4;
u_ind_z = u_z1 + u_z2 + u_z3 + u_z4;
u_ind = [u_ind_x;u_ind_y;u_ind_z];
a_inf = dot(u_ind,vect_norm);
end