%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%       Sparse Pseudo-Spectral Projection
%
%       Author: P. Sochala (version 04-27-20)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


function result = PSP_data(N_VA,level)

n1 = num2str(N_VA) ;
n2 = num2str(level) ;
name1 = ['PC_data/dim',n1,'/level',n2,'/Points.txt'] ;
name2 = ['PC_data/dim',n1,'/level',n2,'/Poids.txt'] ;
name3 = ['PC_data/dim',n1,'/level',n2,'/Indices.txt'] ;
name4 = ['PC_data/dim',n1,'/level',n2,'/Matrix.txt'] ;

%~~~ Points and Weights ~~~~~~~~~~~~~~~~~
PI = load(name1) ;
WI = load(name2) ;
PI = 2*PI-1 ;
N_PI = size(PI,1);

%~~~ Multi-indices of GPC ~~~~~~~~~~~~~~~
Multi_ind = load(name3) ;
N_PC = size(Multi_ind,1) ;

%~~~ Squared norm of GPC ~~~~~~~~~~~~~~~~
c = zeros(1,N_PC) ;
c(1) = 1 ;
for i=2:N_PC
	ind = find(Multi_ind(i,:)) ;
        deg_pol = Multi_ind(i,ind) ;
        c(i) = prod(ones(size(deg_pol,1))./(2*deg_pol+1)) ;
end

%~~~ Matrix Projection ~~~~~~~~~~~~~~~~~~
Coeff = load(name4) ;
A = sparse(Coeff(:,1)+1,Coeff(:,2)+1,Coeff(:,3)) ;
  %      A = A./(repmat(c,N_PI,1).^(0.5))' ;            
for i=1:N_PI
	A(:,i) = A(:,i)./(c.^(0.5))' ;         
end
        
 
result = struct( 'PI'        , PI , ...      
                 'WI'        , WI , ...    
                 'Multi_ind' , Multi_ind , ...
                 'sqnorm'    , c' , ...
                 'Matrix'    , A ) ;


end

