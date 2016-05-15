function [alpha]=blade_angle(span,roottwist,tiptwist,pitch,radial_coord);
    alpha = pitch + roottwist - (roottwist - tiptwist)*(radial_coord/span);
end