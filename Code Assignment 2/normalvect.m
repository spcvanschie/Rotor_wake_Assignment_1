function normalvects = normalvect(alpha,rot_ang)
normalvects = zeros(3,length(alpha));
for i = 1:length(alpha)
    x = sin(deg2rad(alpha(i)));
    y = -ones(1,length(alpha(i))).*sin(rot_ang);
    z = cos(deg2rad(alpha(i)))*cos(rot_ang);
    normal = [x; y; z];
    vect_norm = norm(normal);
    normal = normal/vect_norm;
    normalvects(1,i) = normal(1);
    normalvects(2,i) = normal(2);
    normalvects(3,i) = normal(3);
end
end