% function C_induced = wakeinduction(cp_coords,vect_norm,TE_blade,TE_convected,circulation_history,eps,drag)
% Q_induced = zeros(length(cp_coords(1,:)),length(circulation_history(1,1,:)),length(cp_coords(1,:)));
% for i = 1:length(cp_coords(1,:))
%     for j = 1:length(circulation_history(1,1,:))
%         for k = 1:length(cp_coords(1,:))
%             if j == 1
%                  
%             else
%                 Q_induced(i,j,k) = HSHOE(vect_norm(:,i),cp_coords(:,i),TE_convected(,,length(circulation_history(1,1,:)-(j-1)))
%             end
%         end
%     end
% end
% 
% 
% end




function [Q_induced, Q_ind] = wakeinduction(Q_ind,cp_coords,vect_norm,TE_convected,TE_blade,circulation_history,eps,drag)
Q_induced = zeros(length(cp_coords(1,:)),1);
for i = 1:length(cp_coords(1,:))
    for j = 1:(length(circulation_history(1,1,:)))
        for k = 1:length(cp_coords(1,:))
%             if j+1 > length(circulation_history(1,1,:))
%                 Q_ind(k,j) = HSHOE(vect_norm(:,i),cp_coords(:,i),TE_convected(:,k,j),TE_blade(:,k),TE_blade(:,k+1),TE_convected(:,k+1,j),circulation_history(k,1,j),eps,0);
%             else
            if j == length(circulation_history(1,1,:))
                Q_ind(k,j) = HSHOE(vect_norm(:,i),cp_coords(:,i),TE_convected(:,k,j),TE_blade(:,k),TE_blade(:,k+1),TE_convected(:,k+1,j),circulation_history(k,1,j),eps,drag);
            else
               Q_ind(k,j) = HSHOE(vect_norm(:,i),cp_coords(:,i),TE_convected(:,k,j),TE_convected(:,k,j+1),TE_convected(:,k+1,j+1),TE_convected(:,k+1,j),circulation_history(k,1,j),eps,drag);
            end
        end
    end
    Q_induced(i,1) = sum(sum(Q_ind));
end
end