function z = HT_y(y, fft_BT, downsampling_scale, image_size, s0)
%HT_Y Adjoint operator of H_z:
%     zero-insertion upsampling, followed by transpose blur.
%
% Input:
%   y                  : LR-HSI matrix, size K x (Mh*Nh)
%   fft_BT             : conjugate FFT of spatial blur kernel
%   downsampling_scale : spatial downsampling factor
%   image_size         : [M, N]
%   s0                 : starting index for downsampling
%
% Output:
%   z                  : adjoint result, size K x (M*N)

    M = image_size(1);
    N = image_size(2);

    if nargin < 5
        s0 = 1;
    end

    row_idx = s0:downsampling_scale:M;
    col_idx = s0:downsampling_scale:N;

    Mh = length(row_idx);
    Nh = length(col_idx);

    transpose_flag = false;
    if size(y, 2) == Mh * Nh
        Y = y;                 % K x (Mh*Nh)
    elseif size(y, 1) == Mh * Nh
        Y = y.';               % convert to K x (Mh*Nh)
        transpose_flag = true;
    else
        error('The size of y is inconsistent with the downsampled image size.');
    end

    K = size(Y, 1);
    Z = zeros(K, M * N);

    for k = 1:K
        % Recover the k-th LR spectral band
        band_lr = reshape(Y(k, :), Mh, Nh);

        % Zero-insertion upsampling
        band_up = zeros(M, N);
        band_up(row_idx, col_idx) = band_lr;

        % Adjoint blur operator
        band_adj = real(ifft2(fft2(band_up) .* fft_BT));

        % Convert back to matrix form
        Z(k, :) = band_adj(:).';
    end

    if transpose_flag
        z = Z.';
    else
        z = Z;
    end
end