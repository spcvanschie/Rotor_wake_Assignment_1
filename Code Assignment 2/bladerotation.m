function coords_new = bladerotation(rot_angle,spanwise,coords)
% This is a function to calculate the translation of any point on the rotor
% blade as function of time, depending on the rotation angle
rot_vect = [1; cos(rot_angle); sin(rot_angle)];
coords_new = zeros(3,length(coords(1,:)));
coords_update = rot_vect*spanwise;
for i = 1:length(coords(1,:))
    coords_new(1,i) = coords(1,i);
    coords_new(2,i) = coords_update(2,i);
    coords_new(3,i) = coords_update(3,i);
end
end