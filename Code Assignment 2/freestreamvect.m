function Q_inf = freestreamvect(normalvect,u_inf,u_rotor)
u_combined = u_rotor + u_inf; % combined freestream and local rotor rotation velocity vector
Q_inf = zeros(length(normalvect(1,:)),1);
for i = 1:length(normalvect(1,:))
    Q_inf(i,1) = dot(u_combined(:,i),normalvect(:,i));
end
end