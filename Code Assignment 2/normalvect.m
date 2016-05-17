function normal = normalvect(alpha,rot_ang)
    x = alpha;
    y = -ones(1,length(alpha)).*sin(rot_ang);
    z = cos(deg2rad(alpha))*cos(rot_ang);
    normal = [x; y; z];
end