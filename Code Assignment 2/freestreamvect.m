function Q_inf = freestreamvect(normalvect,u_inf)
Q_inf = zeros(length(normalvect(1,:)),1);
for i = 1:length(normalvect(1,:))
    Q_inf(i,1) = dot(u_inf,normalvect(:,i));
end
end