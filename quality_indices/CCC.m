function [cc,cc_bands] = CCC(X, S)

% CC computes cross-corellation between original dataset and its estimate
% cc = CC(X,S) returns number between 0 and 1
% 
% INPUT ARGUMENTS:
%   X  : recovered SRI data 
%   S  : original SRI data
% OUTPUT ARGUMENTS:
%   cc: cross-correlation between X and SRI



S3 = reshape( S,[size(S,1)*size(S,2), size(S,3)] ); % mode 3 matricization
X3 = reshape( X,  [size(X,1)*size(X,2), size(X,3)] ); 
cc=0;
for k=1:size(S,3)
    cc_bands(k) = corr(S3(:,k),X3(:,k));
    cc = cc + cc_bands(k);    
end
cc = cc/(size(S,3));
cc_bands = cc_bands';
end

