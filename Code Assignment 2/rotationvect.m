function u_rotor = rotationvect(omega,rot_ang,spanwise,normalvect)
v_rot = spanwise./omega; % rotation velocity magnitude of each collocation point
rot_vect = [0; -sin(rot_ang); cos(rot_ang)]; % rotation angle velocity vector
u_rotor = zeros(length(normalvect(1,:)),1);
rot_velocity = rot_vect*v_rot;
for i = 1:length(normalvect(1,:))
    u_rotor(i,1) = dot(rot_velocity(:,i),normalvect(:,i));
end
end