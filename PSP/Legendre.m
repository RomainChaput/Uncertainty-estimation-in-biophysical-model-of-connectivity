%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%     Polynomes de Legendre 1D          
%   
%     Input  : point d'évaluation (X) et degré le plus élevé (ordre)
%     Output : Valeurs des polynômes (rangés par degré croissant) en X
%    
%     Auteur : P.Sochala (dernière version le 19-06-12)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function PL = Legendre(X,ordre) 

PL = zeros(size(X,1),ordre+1) ;

PL(:,1) = ones(size(X,1),1) ;
PL(:,2) = X ;
for i=2:ordre
    PL(:,i+1) = ((2*i-1)*X.*PL(:,i)-(i-1).*PL(:,i-1))/i ;
end

