function [p] = multiLR_slice(MSI_all_lr, H_lr, lbd1, max_iter)
%MULTILR_SLICE Estimate the linear combination coefficients.
%
% Inputs:
%   MSI_all_lr : Downsampled multispectral image, with size [Ih, Jh, K]
%   H_lr       : A low-resolution hyperspectral band, with size [Ih, Jh]
%   lbd1       : Regularization parameter
%   max_iter   : Maximum number of iterations
%
% Output:
%   p          : Estimated coefficient vector, with size [K, 1]

[Ih, Jh, K] = size(MSI_all_lr);
N = Ih * Jh;  % Total number of spatial pixels

%% 1. Convert the 3-D tensor into a 2-D matrix
% Each column corresponds to one spatial pixel, and each row corresponds
% to one multispectral band.
X = reshape(MSI_all_lr, N, K)';  % K-by-N matrix
Y = H_lr(:)';                    % 1-by-N vector

%% 2. Compute the Gram matrix
Gam = X * X';                    % K-by-K matrix

%% 3. Initialize the variables
p = ones(K, 1);                  % Initial coefficient vector
I_reg = lbd1 * eye(K);           % Regularization matrix

%% 4. Iterative optimization
for iter = 1:max_iter

    % Store the coefficient vector from the previous iteration
    p_prev = p;

    % Update the coefficient vector p
    p = (Gam + I_reg) \ (X * Y');

    % Compute the reconstruction error
    H_rec = p' * X;              % 1-by-N reconstructed result
    Res = H_rec - Y;             % 1-by-N residual
    res_norm = norm(Res, 'fro')^2 / N;

    % Display the iteration information
    fprintf('Iteration = %3d, Residual = %.7f\n', iter, res_norm);

    % Check the stopping criterion
    p_change = norm(p - p_prev) / max(norm(p_prev), eps);

    if p_change < 1e-5
        fprintf('The algorithm converged early at iteration %d.\n', iter);
        break;
    end
end

end



