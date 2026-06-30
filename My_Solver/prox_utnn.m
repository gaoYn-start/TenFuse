function [X, tnn, trank] = prox_utnn(UU, Y, rho)
%PROX_UTNN Proximal operator of the transformed tensor nuclear norm.
%
%   Solve:
%       min_X  rho * ||X||_TTNN + 0.5 * ||X - Y||_F^2
%
%   Inputs:
%       UU   - unitary transform matrix along the 3rd mode, size n3 x n3
%       Y    - input tensor, size n1 x n2 x n3
%       rho  - threshold parameter
%
%   Outputs:
%       X     - proximal solution tensor
%       tnn   - transformed tensor nuclear norm of X
%       trank - transformed tubal rank of X
%
%   This implementation does NOT require Tensor Toolbox.

    %% Basic size information
    [n1, n2, n3] = size(Y);

    if size(UU,1) ~= n3 || size(UU,2) ~= n3
        error('The size of UU must be n3-by-n3, where n3 = size(Y,3).');
    end

    %% ---------------------------------------------------------------
    %  Step 1: Mode-3 unfolding of Y
    %  Y_mat has size n3 x (n1*n2)
    %% ---------------------------------------------------------------
    Y_mat = mode3_unfold(Y);

    %% ---------------------------------------------------------------
    %  Step 2: Apply unitary transform along the 3rd mode
    %% ---------------------------------------------------------------
    Y_hat_mat = UU' * Y_mat;

    %% ---------------------------------------------------------------
    %  Step 3: Fold back to tensor form
    %  Y_hat has size n1 x n2 x n3
    %% ---------------------------------------------------------------
    Y_hat = mode3_fold(Y_hat_mat, n1, n2, n3);

    %% ---------------------------------------------------------------
    %  Step 4: Singular value thresholding slice by slice
    %% ---------------------------------------------------------------
    X_hat = zeros(n1, n2, n3);
    tnn   = 0;
    trank = 0;

    for i = 1:n3
        [U, S, V] = svd(Y_hat(:,:,i), 'econ');

        s = diag(S);
        s = max(s - rho, 0);

        r = nnz(s > 0);

        if r > 0
            X_hat(:,:,i) = U(:,1:r) * diag(s(1:r)) * V(:,1:r)';
        else
            X_hat(:,:,i) = zeros(n1, n2);
        end

        tnn   = tnn + sum(s);
        trank = max(trank, r);
    end

    %% ---------------------------------------------------------------
    %  Step 5: Inverse unitary transform along the 3rd mode
    %% ---------------------------------------------------------------
    X_hat_mat = mode3_unfold(X_hat);
    X_mat     = UU * X_hat_mat;

    X = mode3_fold(X_mat, n1, n2, n3);

    %% ---------------------------------------------------------------
    %  Optional normalization
    %  Retained to stay consistent with your original code.
    %% ---------------------------------------------------------------
    tnn = tnn / n3;

end


%% ===================================================================
%  Local function: mode-3 unfolding
%  Input : Tensor A of size n1 x n2 x n3
%  Output: Matrix A_mat of size n3 x (n1*n2)
%% ===================================================================
function A_mat = mode3_unfold(A)
    A_mat = reshape(permute(A, [3, 1, 2]), size(A,3), []);
end


%% ===================================================================
%  Local function: mode-3 folding
%  Input : Matrix A_mat of size n3 x (n1*n2)
%  Output: Tensor A of size n1 x n2 x n3
%% ===================================================================
function A = mode3_fold(A_mat, n1, n2, n3)
    A = permute(reshape(A_mat, [n3, n1, n2]), [2, 3, 1]);
end