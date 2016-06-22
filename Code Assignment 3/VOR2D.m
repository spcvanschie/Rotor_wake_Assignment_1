function vel_ind = VOR2D(gamma,x,z,x_j,z_j,eps)
r_j = (x - x_j)^2 + (z - z_j)^2;
if r_j > eps
    vel_ind = (gamma/(2*pi*(r_j^2)))*[0 1; -1 0]*[x - x_j; z - z_j];
else
    vel_ind = [0;0];
end
end