function P = hermitepols(nord,x)
%
% Returns the Hermite, H_e, polynomials of degree 0 to nord evaluated at the points x
% x must be a column vector
% The output matrix P is of size length(x) x (nord+1) where
% P(i,k+1) = P_k(x_i)
% The polynomials are computed by recurrence.
%
P = ones([length(x) nord+1]);
P(:,1) = 1;
if (nord > 0) 
  P(:,2) = x;
end
for k=2:nord
  P(:,k+1) = x.*P(:,k) - (k-1) * P(:,k-1) ;
end

