function s_bands = PSNR_bands(A, B)
    row = 0;
    col = 0;
    A = double(im2uint8(A));
    B = double(im2uint8(B));
    [n, m, ch] = size(A);
    
    if ch == 1
        % 单通道图像
        e = A - B;
        e_cropped = e(row+1:n-row, col+1:m-col);
        mse = mean(e_cropped(:).^2);  % 使用向量化计算MSE
        s_bands = 10 * log10(255^2 / mse);  % 输出单元素向量
    else
        % 多通道图像：预分配输出向量
        s_bands = zeros(ch,1);
        
        for i = 1:ch
            % 提取当前通道并计算差值
            e = A(:, :, i) - B(:, :, i);
            e_cropped = e(row+1:n-row, col+1:m-col);
            mse = mean(e_cropped(:).^2);  % 使用向量化计算MSE
            s_bands(i) = 10 * log10(255^2 / mse);  % 存储当前通道PSNR
        end
    end
end