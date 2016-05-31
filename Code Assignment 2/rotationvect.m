function rot_velocity = rotationvect(omega,rot_ang,spanwise,normalvect)
v_rot = spanwise./omega; % rotation velocity magnitude of each collocation point
rot_vect = [0; -sin(rot_ang); cos(rot_ang)]; % rotation angle velocity vector
rot_velocity = rot_vect*v_rot;
end