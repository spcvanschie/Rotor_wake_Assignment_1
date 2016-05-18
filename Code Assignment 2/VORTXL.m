function [u_ind]=VORTXL(x_cp,x_a,x_b,circ,eps);
% x_cp, x_a and x_b are all 3D position vectors: 
% x_cp is the position of a random point in space on which the influence of a vortex line element is computed, 
% x_a denotes one end of the vortex line element and x_b the other end
% eps denotes the vortex core size, typically a very small value


% Calculate vortex element vectors
r0_vect = x_b - x_a;
r1_vect = x_cp - x_a;
r2_vect = x_cp - x_b;

r1_cross_r2 = cross(r1_vect,r2_vect);
r1 = norm(r1_vect);
r2 = norm(r2_vect);

if r1<=eps || r2<=eps || norm(r1_cross_r2)^2<=eps
   u_ind = [0; 0; 0]
   return
end

r0_dot_r1 = dot(r0_vect,r1_vect);
r0_dot_r2 = dot(r0_vect,r1_vect);

K = ((r0_dot_r1/r1) - (r0_dot_r2/r2))*circ/(4*pi*(norm(r1_cross_r2)^2));
u_ind = K*r1_cross_r2;

end