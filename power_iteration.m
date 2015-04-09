function [ pHat, tHat ] = power_iteration( A )
% dummy variable
pHat = 0;

%[U,S,V] = svd(A);

%tHat = sign(U(:, 1));
[V,D]= eigs(A*A'); 
tHat = sign(V(:,1)); 

end

