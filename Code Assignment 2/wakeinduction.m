function [Q_induced, Q_ind] = wakeinduction(Q_ind,cp_coords,vect_norm,TE_convected,TE_blade,circulation_history,eps)
Q_induced = zeros(length(cp_coords(1,:)),1);
for i = 1:length(cp_coords(1,:))
    for j = 1:length(circulation_history(1,1,:))
        for k = 1:length(cp_coords(1,:))
            if j+1 > length(circulation_history(1,1,:))
                Q_ind(k,j) = HSHOE(vect_norm(:,i),cp_coords(:,i),TE_convected(:,k,j),TE_blade(:,k),TE_blade(:,k+1),TE_convected(:,k+1,j),circulation_history(k,1,j),eps,0);
            else
                Q_ind(k,j) = HSHOE(vect_norm(:,i),cp_coords(:,i),TE_convected(:,k,j+1),TE_convected(:,k,j),TE_convected(:,k+1,j),TE_convected(:,k+1,j+1),circulation_history(k,1,j),eps,0);
            end
        end
    end
    Q_induced(i,1) = sum(sum(Q_ind,2));
end
end