function rmse = RMSE(X, S)
%--------------------------------------------------------------------------
% root mean square error (RMSE)
%
% USAGE
%   out = RMSE(X, S)
%
% INPUT
%   X  : recovered SRI data 
%   S  : original SRI data
%
% OUTPUT
%   out : RMSE (scalar)
%
%--------------------------------------------------------------------------
[n1,n2,n3]=size(X);
% aux  = sum(sum((X - S).^2, 1), 2)/(n1*n2);
aux  = sum(sum((X - S).^2, 1), 2);
rmse = sqrt(sum(aux, 3)/(n1*n2*n3));
end