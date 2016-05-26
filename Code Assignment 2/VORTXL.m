function [u_x,u_y,u_z]=VORTXL(x_cp,x_a,x_b,circ,eps);
% x_cp, x_a and x_b are all 3D position vectors: 
% x_cp is the position of a random point in space on which the influence of a vortex line element is computed, 
% x_a denotes one end of the vortex line element and x_b the other end
% eps denotes the vortex core size, typically a very small value

r1_cross_r2_x = (x_cp(2,1)-x_a(2,1))*(x_cp(3,1)-x_b(3,1))-(x_cp(3,1)-x_a(3,1))*(x_cp(2,1)-x_b(2,1));
r1_cross_r2_y = -(x_cp(1,1)-x_a(1,1))*(x_cp(3,1)-x_b(3,1))+(x_cp(3,1)-x_a(3,1))*(x_cp(1,1)-x_b(1,1));
r1_cross_r2_z = (x_cp(1,1)-x_a(1,1))*(x_cp(2,1)-x_b(2,1))-(x_cp(2,1)-x_a(2,1))*(x_cp(1,1)-x_b(1,1));

r1_cross_r2_sq = r1_cross_r2_x^2 + r1_cross_r2_y^2 + r1_cross_r2_z^2;
r1 = sqrt((x_cp(1,1)-x_a(1,1))^2+(x_cp(2,1)-x_a(2,1))^2+(x_cp(3,1)-x_a(3,1))^2);
r2 = sqrt((x_cp(1,1)-x_a(1,1))^2+(x_cp(2,1)-x_a(2,1))^2+(x_cp(3,1)-x_a(3,1))^2);

r0_dot_r1 = (x_b(1,1)-x_a(1,1))*(x_cp(1,1)-x_a(1,1))+(x_b(2,1)-x_a(2,1))*(x_cp(2,1)-x_a(2,1))+(x_b(3,1)-x_a(3,1))*(x_cp(3,1)-x_a(3,1));
r0_dot_r2 = (x_b(1,1)-x_a(1,1))*(x_cp(1,1)-x_b(1,1))+(x_b(2,1)-x_a(2,1))*(x_cp(2,1)-x_b(2,1))+(x_b(3,1)-x_a(3,1))*(x_cp(3,1)-x_b(3,1));

K = ((r0_dot_r1/r1)-(r0_dot_r2/r2))*circ/(4*pi*r1_cross_r2_sq);
u_x = K*r1_cross_r2_x;
u_y = K*r1_cross_r2_y;
u_z = K*r1_cross_r2_z;

u_ind = [u_x; u_y; u_z];

if r1<=eps || r2<= eps || r1_cross_r2_sq<= eps
    u_ind = [0;0;0];
end

%% Old approach
% % Calculate vortex element vectors
% r0_vect = x_b - x_a;
% r1_vect = x_cp - x_a;
% r2_vect = x_cp - x_b;
% 
% r1_cross_r2 = cross(r1_vect,r2_vect);
% r1 = norm(r1_vect);
% r2 = norm(r2_vect);
% 
% if r1<=eps || r2<=eps || norm(r1_cross_r2)^2<=eps
%    u_ind = [0; 0; 0]
%    return
% end
% 
% r0_dot_r1 = dot(r0_vect,r1_vect);
% r0_dot_r2 = dot(r0_vect,r1_vect);
% 
% K = ((r0_dot_r1/r1) - (r0_dot_r2/r2))*circ/(4*pi*(norm(r1_cross_r2)^2));
% u_ind = K*r1_cross_r2;

end