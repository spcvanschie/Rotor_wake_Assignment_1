function normal = normalvect(alpha,rot_ang)
    x = alpha;
    y = -ones(1,length(alpha)).*sin(deg2rad(rot_ang));
    z = cos(deg2rad(alpha))*cos(deg2rad(rot_ang));
    normal = [x; y; z];
end