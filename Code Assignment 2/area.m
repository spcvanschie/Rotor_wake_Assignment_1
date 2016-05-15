function A = area(chord,vort_end_y);
A = zeros(1,length(chord));
for i = 1:length(chord);
    A(i) = chord(i)*((vort_end_y(i+1)-vort_end_y(i));
end
end