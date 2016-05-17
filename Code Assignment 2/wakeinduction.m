function Q_ind = wakeinduction(cp_coords,vect_norm,TE_convected,TE_blade,circulation_history,eps)
Q_ind = zeros(length(cp_coords(1,:)),1);

for i = 1:length(cp_coords(1,:))
    for j = 1:length(circulation_history(1,1,:))
        for k = 1:length(cp_coords(1,:))
            if j+1 > length(circulation_history(1,1,:))
                Q_ind(i,1) = Q_ind(i,1) + HSHOE(vect_norm(:,i),cp_coords(:,i),TE_blade(:,k),TE_blade(:,k+1),TE_convected(:,k+1,j),TE_convected(:,k,j),circulation_history(1,k,j),eps);

            else
                Q_ind(i,1) = Q_ind(i,1) + HSHOE(vect_norm(:,i),cp_coords(:,i),TE_convected(:,k,j+1),TE_convected(:,k+1,j+1),TE_convected(:,k+1,j),TE_convected(:,k,j),circulation_history(1,k,j),eps);
            end
        end
    end
end
end