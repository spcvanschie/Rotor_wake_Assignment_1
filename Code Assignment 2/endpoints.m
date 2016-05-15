function [points] = endpoints(root,tip,cp_all);
points = zeros(1,length(cp_all)+1);
points(1) = root;
points(length(cp_all)+1) = tip;

for i = 2:(length(cp_all));
    points(i) = (cp_all(i-1)+cp_all(i))/2;
end
end