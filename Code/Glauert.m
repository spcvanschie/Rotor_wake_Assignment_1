function C_t=Glauert(a)
C_t = ones(1,length(a));
C_t1 = 1.816;
a_1 = 1 + 0.5*sqrt(C_t1);
for i = (1:length(a))
    if a(i) >= a_1
        C_t(i) = C_t1 - 4*(sqrt(C_t1)-1)*(1-a(i));
    else
        C_t(i) = 4*a(i)*(1-a(i));
    end
end
end