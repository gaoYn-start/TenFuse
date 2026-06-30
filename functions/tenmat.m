function O = tenmat(Y, row_modes)
%TENMAT Lightweight tensor matricization function.
%
%   O = tenmat(Y, [n])
%
%   This simplified implementation supports tensor unfolding with
%   specified row modes. For a 3-D tensor Y, tenmat(Y,[3]) returns
%   the mode-3 unfolding matrix of size:
%
%       size(Y,3) x (size(Y,1)*size(Y,2))
%
%   This version is sufficient for the use in prox_utnn.m.

    dims = size(Y);
    N = ndims(Y);

    if nargin < 2
        error('tenmat requires at least two inputs.');
    end

    row_modes = row_modes(:).';
    col_modes = setdiff(1:N, row_modes, 'stable');

    % Reorder tensor dimensions: row modes first, column modes second
    Y_perm = permute(Y, [row_modes, col_modes]);

    % Compute matrix size
    row_size = prod(dims(row_modes));
    col_size = prod(dims(col_modes));

    % Matricization
    O = reshape(Y_perm, row_size, col_size);
end