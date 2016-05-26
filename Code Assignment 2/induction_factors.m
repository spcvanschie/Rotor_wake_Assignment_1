function A_ind = induction_factors(vect_norm,x_cp,vort_end,TE,eps,drag);
A_ind = zeros(length(x_cp(1,:)),length(x_cp(1,:)));
for i = 1:length(x_cp(1,:))
    for j = 1:length(x_cp(1,:))
        A_ind(i,j) = HSHOE(vect_norm(:,i),x_cp(:,i),vort_end(:,j),vort_end(:,j+1),TE(:,j),TE(:,j+1),1,eps,drag);
    end
end