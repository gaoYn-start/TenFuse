function [rgb] = func_hyperImshowR_2( hsi, RGBbands )
%% Hyperspectral Image color display
% Author: Zephyr Hou
% Time: 2019-12-02
% Modified to enhance brightness

hsi = double(hsi);
[rows, cols, bands] = size(hsi);

% 使用分位数截断归一化以避免离群值影响
minVal = prctile(hsi(:), 2); % 2%分位数作为最小值
maxVal = prctile(hsi(:), 98); % 98%分位数作为最大值
normalizedData = (hsi - minVal) / (maxVal - minVal);
normalizedData = max(0, min(normalizedData, 1)); % 截断到[0,1]范围

hsi = normalizedData;

if (nargin == 1)
    RGBbands = [bands round(bands/2) 1];
end

% 提取RGB波段
if (bands == 1)
    red = hsi(:,:);
    green = hsi(:,:);
    blue = hsi(:,:);
else
    red = hsi(:,:,RGBbands(1));
    green = hsi(:,:,RGBbands(2));
    blue = hsi(:,:,RGBbands(3));
end

% 调整自适应直方图均衡化参数并增加亮度增益
clipLimit = 0.005;    % 增大对比度限制
gammaValue = 0.95;     % Gamma校正参数（小于1提升亮度）

rgb = zeros(rows, cols, 3);
rgb(:,:,1) = adapthisteq(red, 'ClipLimit', clipLimit).^gammaValue;
rgb(:,:,2) = adapthisteq(green, 'ClipLimit', clipLimit).^gammaValue;
rgb(:,:,3) = adapthisteq(blue, 'ClipLimit', clipLimit).^gammaValue;

% 可选：进一步线性增强亮度
brightnessGain = 1.2; % 亮度增益系数
rgb = min(rgb * brightnessGain, 1); % 限制最大值不超过1

imshow(rgb); axis image;
end