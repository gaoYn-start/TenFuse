function[d]=row_extraction(b31k) 
n_rows = length(b31k);

% 生成n_rows阶单位矩阵
I = eye(n_rows);

% 确定要选取的行索引
if mod(n_rows, 2) == 0
    idx = n_rows / 2;      % 偶数行数，取中间的上一行
else
    idx = (n_rows + 1) / 2; % 奇数行数，取正中间一行
end

% 提取对应的行向量
d = I(idx, :);