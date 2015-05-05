function Za_hat = normalizeMatrix(Zaj, Zbj, Zcj)

% make sure that all labels average vecotrs have the same length (i.e. number of items)
assert(length(Zaj) == length(Zbj)); 
assert(length(Zcj) == length(Zbj)); 

n = length(Zaj); 

A = zeros(size(Zaj,2)); 
B = zeros(size(Zaj,2));
for i = 1:n 
    A = Zaj(i,:)' * Zbj(i,:); 
    B = Zbj(i,:)' * Zcj(i,:); 
end
A = A / n; 
B = B / n; 

Za_hat = zeros(n,size(Zaj,2));
tmp_variable = B * (A^(-1)); 
for i = 1:n
    Za_hat(i,:) = tmp_variable * Zaj(i,:)'; 
end 
