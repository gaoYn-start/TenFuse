function[M] = tenv3(MSI,p)


[I,J,K] = size(MSI);
M       = zeros(I,J);

for i = 1:K
    deM = MSI(:,:,i)*p(i);
    M   = M +deM;
    
end
