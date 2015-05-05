clear all; close all; clc; 

m = 100; 
n = 200; 
%Z = unidrnd(2,m,n) -1;   % random matrix of 0 and 1 
Z = zeros(m, n, 2); 
for i = 1:m 
    for j = 1:n
        ind = unidrnd(2); 
        Z(i, j, ind) = 1; 
    end
end

% partition the Z matrix into three parts 
inds = unidrnd(3, m, 1); 
Za = Z(inds==1,:,:); 
Zb = Z(inds==2,:,:); 
Zc = Z(inds==3,:,:); 

% average label of each item 
Zaj_tmp = sum(Za)/size(Za,1); 
Zbj_tmp = sum(Zb)/size(Zb,1); 
Zcj_tmp = sum(Zc)/size(Zc,1); 
Zaj = reshape(Zaj_tmp, size(Zaj_tmp,2), size(Zaj_tmp,3)); 
Zbj = reshape(Zbj_tmp, size(Zbj_tmp,2), size(Zbj_tmp,3)); 
Zcj = reshape(Zcj_tmp, size(Zcj_tmp,2), size(Zcj_tmp,3)); 

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

AA = rand(10,1);
[U,S,V] = svd(AA*AA'); 

schur(AA*AA')

AA = rand(10,1);
AAA = AA*AA';  % random positive definte matrix 
[U S V]=svd(AAA);
d=diag(U'*AAA*U);
S=S*spdiags(sign(d.*diag(S)),0,size(U,2),size(U,2));
UU = U * (diag(diag(S)))^(0.5); 

Areconstruct1 = U*S*U'; 
Areconstruct2 = UU*eye(size(U,2))*UU'; 

% Check
norm(Areconstruct1-Areconstruct2)
norm(Areconstruct1-AAA)/norm(AAA)
norm(Areconstruct2-AAA)/norm(AAA)

Q = UU^(-1); 
Q*AAA*Q'; 


[U,S,V] = svd(M2_hat); 
Q = U;

% calculate the whitened tensor
for i = 1:size(M3_hat,1)
    for j = 1:size(M3_hat,1)
        for k = 1:size(M3_hat,1)
            M_3QQQ(i,j,k) = 0; 
            for l = 1:size(M3_hat,1)
                for m = 1:size(M3_hat,1)
                    for t = 1:size(M3_hat,1)
                        M_3QQQ(i,j,k) = M_3QQQ(i,j,k) + ... 
                            M3_hat(l,m,t) * Q(l,i) * Q(m,j) * Q(t,k); 
                    end
                end
            end
        end
    end
end


%K - Number of hidden components
params.K = 2; 
%E - Number of power iterations
params.E = 10;  
%v - Variance of the observation model, see the tensor factorization paper for details
params.v = 1; 
%display - Toggle to display iterations, if 1 code shows the iterations
params.display = 1; 

[muhat, what] = tensor_power_iterations(M_3QQQ,params); 

[muhat_maxVal, muhat_maxInd] = max(muhat); 
C_hat(muhat_maxInd) = muhat;
