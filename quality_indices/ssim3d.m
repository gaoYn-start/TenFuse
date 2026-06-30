function [ssim]  =  ssim3d(imagery1,imagery2)
ssim = 0;
for i = 1:size(imagery1,3)
    ssim = ssim + ssim_index(imagery1(:, :, i),imagery2(:, :, i));
end
ssim= ssim/size(imagery1,3);
end

