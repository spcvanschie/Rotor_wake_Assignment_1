function [a_inf] = HSHOE(vect_norm,x_cp,x_a,x_b,x_c,x_d,circ,eps,drag);
u_1 = VORTXL(x_cp,x_a,x_b,circ,eps);
u_2 = VORTXL(x_cp,x_b,x_c,circ,eps);
u_3 = VORTXL(x_cp,x_c,x_d,circ,eps);
u_4 = VORTXL(x_cp,x_d,x_a,circ,eps);
if drag > 0:
    u_2 = [0;0;0];
    u_4 = [0;0;0];
end
u_ind = u_1 + u_2 + u_3 + u_4;
a_inf = dot(u_ind,vect_norm);
end