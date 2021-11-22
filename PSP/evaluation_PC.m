
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     Evaluation des Polynômes de Chaos MultiD en des points        
%   
%     Input  : Points et Multi-indices des Polynômes
%     Output : GPC aux points 
%    
%     Auteur : P.Sochala (dernière version le 18-12-14)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function GPC = evaluation_PC(Points , Multi_ind)

N_VA = size(Points,2) ;                     % Number of RV
ordre = max(Multi_ind(:)) ;                 % Maximal order 

%~~~ Old version ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% N_PC = size(Multi_ind,1) ;                  % Number of GPC
% PL_ND = zeros(size(Points,1),N_PC,N_VA) ;
% 
% for i=1:N_VA
%     PL_1D = Legendre(Points(:,i),ordre) ;
%     PL_ND(:,:,i)=PL_1D(:,Multi_ind(:,i)+1) ;
% end
% 
% GPC = prod(PL_ND,3) ;

%~~~ New version ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PL_1D = Legendre(Points(:,1),ordre) ;
GPC=PL_1D(:,Multi_ind(:,1)+1) ;

% size(GPC)
% size(PL_1D)
%return

for i=2:N_VA
    PL_1D = Legendre(Points(:,i),ordre) ;
    indices = find(Multi_ind(:,i)) ;
    GPC(:,indices) = GPC(:,indices).*PL_1D(:,Multi_ind(indices,i)+1) ;
%~~~ Version "dévectorisée" pour les niveaux élevés ~~~~~~~~~~~~~~~~~~~~ 
%     for j=1:size(indices,1)
%         GPC(:,indices(j)) = GPC(:,indices(j)).*PL_1D(:,Multi_ind(indices(j),i)+1) ;
%     end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end

