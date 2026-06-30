function S_sharpened = InitialS_band(MSI, HSI, ratio, lbd1)
%% Input: 
%%% MSI: High-resolution multispectral image [I, J, Km]
%%% HSI: Low-resolution hyperspectral image [Ih, Jh, Kh]
%%% ratio: Resolution ratio (e.g., 4 means that the spatial resolution
%%%        of HSI is one-fourth of that of MSI)
%%% lbd1: Regularization parameter

[I, J, Km] = size(MSI);
[Ih, Jh, Kh] = size(HSI);

% Verify the resolution ratio
if Ih ~= I/ratio || Jh ~= J/ratio
    error('Incorrect resolution ratio: I/ratio=%d but Ih=%d, J/ratio=%d but Jh=%d',...
          I/ratio, Ih, J/ratio, Jh);
end

% 1. Prepare the data
de_Msi = imresize(MSI, [Ih, Jh], 'bicubic');  % Downsample to the low resolution
de_Msi_all = cat(3, de_Msi, ones(Ih, Jh));    % Add a constant term [Ih, Jh, Km+1]
MSI_all = cat(3, MSI, ones(I, J));             % High-resolution data [I, J, Km+1]

% 2. Preallocate the output
S_sharpened = zeros(I, J, Kh);
max_iter = 10;

% 3. Process each HSI band
t1 = tic;
for i = 1:Kh
    H_band = HSI(:, :, i);  % Current low-resolution band
    
    % Call the regression function
    [p] = multiLR_slice(de_Msi_all, H_band, lbd1, max_iter);
    
    % Reconstruct the high-resolution band using the coefficients
    band_rec = zeros(I, J);
    for c = 1:Km+1
        band_rec = band_rec + p(c) * MSI_all(:, :, c);
    end
    S_sharpened(:, :, i) = band_rec;
end

t2 = toc(t1);
fprintf('ALM took %2.2f seconds to run.\n', t2);
end























