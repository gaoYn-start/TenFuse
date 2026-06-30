function  result  = imblizoom( jpg,Multiple )
    %%
	[row,clounm,ID] = size(jpg);
	row1 = round(row*Multiple); % 计算缩放后的图像高度，最近取整
	clounm1 = round(clounm*Multiple); % 计算缩放后的图像宽度，最近取整
	result = zeros(row1,clounm1,ID); % 创建新图像
	%%  扩展矩阵
	expand = zeros(row+2,clounm+2,ID);
	expand(2:row+1,2:clounm+1,:) = jpg;
	expand(1,2:clounm+1,:)=jpg(1,:,:);expand(row+2,2:clounm+1,:)=jpg(row,:,:);
	expand(2:row+1,1,:)=jpg(:,1,:);expand(2:row+1,clounm+2,:)=jpg(:,clounm,:);
	expand(1,1,:) = jpg(1,1,:);expand(1,clounm+2,:) = jpg(1,clounm,:);
	expand(row+2,1,:) = jpg(row,1,:);
	expand(row+2,clounm+2,:) = jpg(row,clounm,:);
	%%  
	for zj = 1:clounm1         % 对图像进行按列逐元素扫描
	    for zi = 1:row1
	        ii = (zi-1)/Multiple; jj = (zj-1)/Multiple;
	        i = floor(ii); j = floor(jj); % 向下取整
	        u = ii - i; v = jj - j;
	        i = i + 1; j = j + 1;
	        result(zi,zj,:) = (1-u)*(1-v)*expand(i,j,:) +(1-u)*v*expand(i,j+1,:) + u*(1-v)*expand(i+1,j,:) +u*v*expand(i+1,j+1,:);
	    end
	end
%	result = uint8(result);
    
end  