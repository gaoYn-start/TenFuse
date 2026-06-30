function rsnr = RSNR(X, S)
%--------------------------------------------------------------------------
% reconstruction signal to noise ratio (RSNR)
%
% USAGE
%   out = RSNR(X, S)
%
% INPUT
%   X  : recovered SRI data 
%   S  : original SRI data
%
% OUTPUT
%   out : RSNR (scalar)
%
%--------------------------------------------------------------------------

num = 0; denum = 0;
for k = 1:size(S,3)
    num   = num + norm(S(:,:,k),'fro')^2;
    denum = denum + norm(S(:,:,k)-X(:,:,k),'fro')^2;
end
rsnr = 10*log10(num/denum);
end