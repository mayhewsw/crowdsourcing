function [ pHat, tHat, T ] = simplified_bp( A )
%Simplified BP
%   As described in class.

% a dummy variable.
pHat = 0;

% max number of iterations
T = 10;

[n,m] = size(A);

% x,y is initialized for every edge i,j
x = zeros(n,m);
y = normrnd(1,1,m,n);

% actually do the iteration updates
for iter=1:T
    
    % update x
    for i=1:n
        for j=1:m
            x(i,j) = A(i, :) * y(:, i) - A(i,j)*y(j,i);
        end
    end
    
    % update y
    for i=1:n
        for j=1:m
            y(j,i) = A(:, j)' * x(:, j) - A(i,j)*x(i,j);
        end
    end
    
end

tHat = zeros(n,1);
for i=1:n
    tHat(i) = sign(A(i, :) * y(:,i));
end


end

