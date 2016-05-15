function A = surfarea(chordlength,vort_y)
A = zeros(1,length(chordlength));
for i = 1:length(chordlength);
    A(i) = chordlength(i)*((vort_y(i+1)-vort_y(i)));
end
end