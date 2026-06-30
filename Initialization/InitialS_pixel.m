function S_sharpened = InitialS_pixel(MSI, HSI, ratio, lbd1)
% Input: 
%   MSI: High-resolution multispectral image [I, J, Km]
%   HSI: Low-resolution hyperspectral image [Ih, Jh, K]
%   ratio: Resolution ratio (e.g., 4 means that the spatial resolution of
%          HSI is one-fourth of that of MSI)
%   lbd1, lbd2: Regularization parameters

[I, J, Km] = size(MSI);
[Ih, Jh, K] = size(HSI);

% Verify the compatibility of the input dimensions
if Ih ~= I/ratio || Jh ~= J/ratio
    error('The resolution ratio is incompatible with the input image dimensions');
end

% Prepare the data matrices with a constant term
de_Msi = imresize(MSI, [Ih, Jh], 'bicubic');  % Downsample MSI to the HSI resolution
de_Msi_all = cat(3, de_Msi, ones(Ih, Jh));    % Add a constant term [Ih, Jh, Km+1]
MSI_all = cat(3, MSI, ones(I, J));             % Add a constant term [I, J, Km+1]

% Preallocate the coefficient tensor
Constant_P = zeros(K, Ih, Jh, Km+1);

% Compute the regression coefficients on the low-resolution grid
for k = 1:K
    H_band = HSI(:, :, k);  % Current low-resolution hyperspectral band
    
    for ih = 1:Ih
        for jh = 1:Jh
            % Obtain the current low-resolution pixel
            h = H_band(ih, jh);
            m_vec = squeeze(de_Msi_all(ih, jh, :));
            
            % Solve the ridge regression problem
            % Replace the missing multiLR_pixel function
            P = (m_vec * m_vec' + lbd1 * eye(Km+1)) \ (m_vec * h);
            
            % Store the coefficients
            Constant_P(k, ih, jh, :) = P;
        end
    end
end

% ====== Reconstruct the high-resolution HSI ======
HSI_sharpened = zeros(I, J, K);

% Upsample the coefficients and reconstruct the image
for k = 1:K
    % Create temporary storage
    band_reconstructed = zeros(I, J);
    
    % Process each coefficient channel separately
    for c = 1:Km+1
        % Upsample the current coefficient channel
        coeff_band = squeeze(Constant_P(k, :, :, c));
        upsampled_coeff = imresize(coeff_band, [I, J], 'bicubic');
        
        % Apply the contribution of the corresponding MSI channel
        band_reconstructed = band_reconstructed + ...
            upsampled_coeff .* MSI_all(:, :, c);
    end
    
    HSI_sharpened(:, :, k) = band_reconstructed;
end

% Output the result
S_sharpened = HSI_sharpened;

% Optional: Display the processing time
% toc;
end
```
