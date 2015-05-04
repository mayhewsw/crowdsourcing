clear all; close all; clc; 

m = 100; 
n = 200; 
Z = unidrnd(2,m,n) -1;   % random matrix of 0 and 1 

% partition the Z matrix into three parts 
inds = unidrnd(3, m, 1); 
Za = Z(inds==1,:); 
Zb = Z(inds==2,:); 
Zc = Z(inds==3,:); 

% average label of each item 
Zaj = sum(Za)/size(Za,1); 
Zbj = sum(Zb)/size(Zb,1); 
Zcj = sum(Zc)/size(Zc,1); 

Zaj_hat = normalizeMatrix(Zaj, Zbj, Zcj); 
Zbj_hat = normalizeMatrix(Zbj, Zaj, Zcj); 

M2_hat = 0; 
M3_hat = 0; 
for i = 1:n 
    M2_hat = M2_hat + Zaj_hat(i) * Zbj_hat(i);    
    M3_hat = M3_hat + Zaj_hat(i) * Zbj_hat(i) * Zcj(i); 
end
M2_hat = M2_hat / n; 
M3_hat = M3_hat / n; 

[U,S,V] = svd(M2_hat); 



