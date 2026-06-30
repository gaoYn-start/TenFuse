function [s] = norm_sparse(X)


Ss = svds(X, 1); 
s = Ss(1);




end