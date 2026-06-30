function ergas = ERGAS(X,S,ratio)
%--------------------------------------------------------------------------
% Erreur relative globale adimensionnelle de synthese (ERGAS)
%
% USAGE
%   ergas = ERGAS(X,S,ratio)
%
% INPUT
%   X  : recovered SRI data 
%   S  : original SRI data
%   ratio : resize factor
%
% OUTPUT
%   ERGAS : ergas (scalar)
%--------------------------------------------------------------------------

ergas_temp = 0;
for k = 1:size(S,3)
    mu = mean(mean(S(:,:,k)));
    ergas_temp = ergas_temp + norm(S(:,:,k)-X(:,:,k),'fro')^2/mu^2;
end
ergas = ratio*sqrt( ergas_temp/numel(S) );
end