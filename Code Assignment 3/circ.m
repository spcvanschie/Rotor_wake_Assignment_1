function [a,RHS,gamma] = circ(timestep,n,V_inf,cp_coords,vort_coords,normal,alpha,omega,gamma_previous,wake_coords,eps)
a = zeros(n,n);
RHS = zeros(n,1);
for j = 1:n % vortex number
    for i = 1:n % collocation point number
        vor = VOR2D(1,cp_coords(1,i),cp_coords(2,i),vort_coords(1,j),vort_coords(2,j),eps);
        a(i,j) = dot(vor,normal);
    end
end
for i = 1:n % collocation point number
    wake_ind = [0;0];
    for k = 1:timestep
        wake_ind = wake_ind + VOR2D(1,cp_coords(1,i),cp_coords(2,i),wake_coords(1,k),wake_coords(2,k),eps);
    end
    v_ind = [cos(deg2rad(alpha)) -sin(deg2rad(alpha)); sin(deg2rad(alpha)) cos(deg2rad(alpha))]*[V_inf(1);V_inf(2)] + [0; omega*cp_coords(1,i)];
    RHS(i,1) = -dot((v_ind+wake_ind),normal);
end


gamma_shed_previous = sum(gamma_previous) - gamma_previous(end);
a = cat(2,a,zeros(n,1));
a = cat(1,a,ones(1,n+1));
RHS = cat(1,RHS,gamma_shed_previous);

for l = 1:n % collocation point number
    vor = VOR2D(1,cp_coords(1,l),cp_coords(2,l),wake_coords(1,timestep + 1),wake_coords(2,timestep + 1),eps);    
    a(l,n+1) = dot(vor,normal);
end


gamma = RHS\a;

end