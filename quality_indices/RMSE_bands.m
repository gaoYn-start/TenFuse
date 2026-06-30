function rmse_bands = RMSE_bands(X, S)
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
Error = X - S;

for i = 1:n3
    aux(i) = sum(sum(Error(:,:,i).^2, 1), 2);
    rmse_bands(i) = sqrt(aux(i)/(n1*n2));
end
rmse_bands = rmse_bands';