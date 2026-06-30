function [sam, map] = SAM(X, S)
%--------------------------------------------------------------------------
% Spectral angle mapper (SAM)
%
% USAGE
%   [angle_SAM, map] = SAM(X, S)
%
% INPUT
%   X  : recovered SRI data 
%   S  : original SRI data
%
% OUTPUT
%   angle_SAM : sam (scalar)
%   map       : 2-D 
%
%--------------------------------------------------------------------------

prod_scal = dot(X,S,3); 
norm_orig = dot(S,S,3);
norm_fusa = dot(X,X,3);
prod_norm = sqrt(norm_orig.*norm_fusa);
prod_map  = prod_norm;
prod_map(prod_map==0)=eps;
map = acos(prod_scal./prod_map);
sam = mean(map(:));
end