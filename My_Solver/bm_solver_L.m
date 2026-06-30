function[blur_ori,P_new]=bm_solver_L(N,H,D,blur_k)


[I,~] = size(N);

%% Transformation into Fourier domain:

pk  = psf2otf(blur_k,[I 1]);
%% Computing the gradient:
F = dftmtx(I)/sqrt(I); 
% Hfft = fft(H);
% H = K*F*H;
B     = F'*D'*D*F;
C     = F'*N*N'*F;
E     = F'*D'*H*N'*F;
EE    = B*diag(pk)*C - E;
nabla_d  = diag(EE);
%% Computing the stepsize:
hessian_d = F'*D'*D*N*N'*F;
alpha = 1/norm(hessian_d,2);
%% GD in Fourier domain:
p_new = pk - alpha*nabla_d;

%% In original domain:
blur_ori_new = otf2psf(p_new);


% if mod(I,2)==0
%    n        = (I-length(blur_k)-1)/2;
%    blur_ori = blur_ori_new(n+2:n+1+length(blur_k),:);
% else
%    n        = (I-length(blur_k))/2;
%    blur_ori = blur_ori_new(n+1:n+length(blur_k),:);
% end
blur_ori = blur_ori_new;
blur_ori(blur_ori<0)=0;
blur_ori = blur_ori/sum(blur_ori);


%% Generating new circulant matrix
p_ori = psf2otf(blur_ori,[I 1]);
C_new = F*diag(p_ori)*F';
C_new = real(C_new);
C_new(C_new<1e-10)=0;

%% Generate new matrix
P_new = D*C_new;










