function [ pHat, tHat ] = power_iteration( A )
% dummy variable
pHat = 0;

% [U,S,V] = svd(A);
% tHat = -sign(U(:, 1));

[V2,D]= eigs(A*A');
D(1,1)

tHat = -sign(sign(sum(V2(:,1)))*V2(:,1)); 

end

