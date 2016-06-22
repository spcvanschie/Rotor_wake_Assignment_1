function [panel_TE,vort_chord,cp_chord,vort_coords,cp_coords,TE_coords] = coordinates(n,c,alpha)
panel_TE = linspace(0,c,n+1);
vort_chord = zeros(1,n);
cp_chord = zeros(1,n);
vort_coords = zeros(2,n);
cp_coords = zeros(2,n);
for i = 2:length(panel_TE)
    vort_chord(i-1) = (panel_TE(i) - panel_TE(i-1))*0.25 + panel_TE(i-1);
    cp_chord(i-1) = (panel_TE(i) - panel_TE(i-1))*0.75 + panel_TE(i-1);
    vort_coords(:,i-1) = vort_chord(i-1)*[cos(deg2rad(alpha));-sin(deg2rad(alpha))];
    cp_coords(:,i-1) = cp_chord(i-1)*[cos(deg2rad(alpha));-sin(deg2rad(alpha))];
end
TE_coords = [c*cos(deg2rad(alpha)); -c*sin(deg2rad(alpha))];

end