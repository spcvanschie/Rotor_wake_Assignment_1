function Thrust_convergence=Thrust_convergence_history(Thrust_history,Thrust)
[height,width] = size(Thrust_history);
for j = (1:height)
    for i = (1:width)
       if Thrust_history(j,i) == 0
           Thrust_history(j,i) = Thrust(i);
       end
    end
end
Thrust_convergence = sum(Thrust_history,2);
end