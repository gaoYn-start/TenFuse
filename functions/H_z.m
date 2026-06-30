function y = H_z(z, fft_B, downsampling_scale, image_size, s0)
%H_Z Spatial degradation operator:
%    blur each band by circular convolution, then downsample.
%
% Input:
%   z                  : 2D hyperspectral matrix, size K x (M*N)
%   fft_B              : FFT of spatial blur kernel
%   downsampling_scale : spatial downsampling factor
%   image_size         : [M, N]
%   s0                 : starting index for downsampling
%
% Output:
%   y                  : degraded LR-HSI matrix, size K x (M/s)*(N/s)

    M = image_size(1);
    N = image_size(2);

    if nargin < 5
        s0 = 1;
    end

    % Automatically handle two possible matrix layouts
    transpose_flag = false;
    if size(z, 2) == M * N
        Z = z;                 % K x (M*N)
    elseif size(z, 1) == M * N
        Z = z.';               % convert to K x (M*N)
        transpose_flag = true;
    else
        error('The size of z is inconsistent with image_size = [%d, %d].', M, N);
    end

    K = size(Z, 1);

    row_idx = s0:downsampling_scale:M;
    col_idx = s0:downsampling_scale:N;

    Mh = length(row_idx);
    Nh = length(col_idx);

    Y = zeros(K, Mh * Nh);

    for k = 1:K
        % Recover the k-th spectral band
        band = reshape(Z(k, :), M, N);

        % Circular convolution by FFT
        band_blur = real(ifft2(fft2(band) .* fft_B));

        % Spatial downsampling
        band_lr = band_blur(row_idx, col_idx);

        % Convert back to 2D matrix form
        Y(k, :) = band_lr(:).';
    end

    if transpose_flag
        y = Y.';
    else
        y = Y;
    end
end