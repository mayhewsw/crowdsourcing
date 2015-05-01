function [ pHat, tHat] = hits( A )
%Simplified BP
%   As described in class.

% a dummy variable.
pHat = 0;

[n,m] = size(A);

E = zeros(2*n, m);
% rewrite A, given that it is binary.
for row=1:n
    Arow = A(row, :);
    top = zeros(1,m);
    top(Arow==-1) = 1;
    
    bottom = zeros(1,m);
    bottom(Arow==1) = 1;
    
    E(2*row-1, :) = top;
    E(2*row, :) = bottom;
end

% E is the adjacency matrix, remains unchanged. 
A = E;

n = 2*n;

for iter=1:10
    y = sum(A, 2);
    y = y / sum(y);
    A = E.*repmat(y, 1, m);
    
    x = sum(A);
    x =x / sum(x);
    A = E.*repmat(x, n, 1);
    
end

n = n/2;

tHat = zeros(n,1);

for i=1:n
    top = y(2*i-1);
    bottom = y(2*i);
    
    if top == 0 && bottom == 0
        tHat(i) = 0;
    elseif top > bottom
        tHat(i) = -1;
    elseif top < bottom
        tHat(i) = 1;
    end
end


end
