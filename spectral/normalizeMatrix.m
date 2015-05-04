function Za_hat = normalizeMatrix(Zaj, Zbj, Zcj)

% make sure that all labels average vecotrs have the same length (i.e. number of items)
assert(length(Zaj) == length(Zbj)); 
assert(length(Zcj) == length(Zbj)); 

n = length(Zaj); 

A = zeros(1,n); 
B = zeros(1,n);
for i = 1:n 
    A(i) = Zaj(i) * Zbj(i); 
    B(i) = Zbj(i) * Zcj(i); 
end
A = A / n; 
B = B / n; 

Za_hat = zeros(1,n); 
for i = 1:n
    Za_hat(i) = B(i) * (A(i)^(-1)) * Zaj(i); 
end 
