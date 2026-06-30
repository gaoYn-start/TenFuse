function Out = quality_indices(X,S,ratio)
%--------------------------------------------------------------------------
% Quality Indices
%
% USAGE
%   Out = quality_indices(X,S,ratio)
%
% INPUT
%   X  : recovered SRI data 
%   S  : original SRI data
%   ratio : GSD ratio between HS and MS imagers
%
% OUTPUT
%   Out.rmse : rmse
%   Out.rsnr : rsnr
%   Out.psnr : psnr
%   Out.ssim : ssim
%   Out.sam  : sam
%   Out.cc   : cc
%   Out.ergas: ergas
%   Out.uiqi : uiqi
%--------------------------------------------------------------------------  

    %%%% Rsnr
    Out.rsnr = RSNR(X,S);  
    %%% Psnr
    Out.psnr = PSNR(X,S);
    %%%% SSIM
    Out.ssim = ssim3d(X*255, S*255);
    %%%% CC
    [cc,cc_bands] = CCC(X,S);
    Out.cc = mean(cc);
    Out.cc_bands = cc_bands;
    %%%% UIQI
    q_band   = zeros(1, size(X,3));
    for idx1 = 1:size(X,3)
        q_band(idx1)=img_qi(S(:,:,idx1), X(:,:,idx1));
    end
    Out.uiqi = mean(q_band);
    %%%% Rmse
    Out.rmse = RMSE(X, S);
    Out.rmse_bands = RMSE_bands(X, S);
    %%%% SAM
    [sam,map] = SAM(X,S);
    Out.sam_map = map;
    Out.sam = sam;
    %%%% ERGAS
    Out.ergas = ERGAS(X,S,ratio);
    


end
