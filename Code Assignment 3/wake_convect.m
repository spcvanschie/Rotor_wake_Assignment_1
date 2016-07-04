function [wake_convected,wake_movement,wake_change,wake_convection] = wake_convect(n,timestep,dt,wake_coords,gamma,vortex_strength,vort_coords,V_inf,eps)
wake_movement = zeros(2,timestep+1);
for k = 1:timestep+1
    if k < timestep + 1
        for i = 1:n
            wake_movement(:,k) = wake_movement(:,k) + VOR2D(gamma(i),wake_coords(1,k),wake_coords(2,k),vort_coords(1,i),vort_coords(2,i),eps);
        end
        for j = (1):(timestep+1)
            vor = VOR2D(vortex_strength(j),wake_coords(1,k),wake_coords(2,k),wake_coords(1,j),wake_coords(2,j),eps);
            wake_movement(:,k) = wake_movement(:,k) + vor;
        end
    end
    wake_movement(:,k) = wake_movement(:,k) + V_inf;
end
wake_change = dt.*wake_movement;
wake_convection = cat(2,wake_change,zeros(2,(length(wake_coords(1,:))-(timestep+1))));
wake_convected = wake_coords + wake_convection;
end