function normal = normalvect(alpha)
n = [sin(deg2rad(alpha));cos(deg2rad(alpha))];
n_x = n(1)/norm(n);
n_z = n(2)/norm(n);
normal = [n_x,n_z];
end