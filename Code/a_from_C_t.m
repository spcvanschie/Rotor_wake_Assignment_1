function a = a_from_C_t(C_t)
    C_t1 = 1.816;
    a_trans = 1 - 0.5*sqrt(C_t1);
    a_linear = (C_t1 - C_t +4 - 4*sqrt(C_t1))/(4-4*sqrt(C_t1)); % Glauert-corrected value for a
    a_BEM = min(roots([-4 4 -C_t])); % BEM-predicted value for a
    if a_BEM > a_trans
        a = a_linear;
    else
        a = a_BEM;
    end
end