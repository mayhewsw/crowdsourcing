function [ pHat, tHat ] = power_iteration( A )
% dummy variable
pHat = 0;

[V,D] = eigs(A');

tHat = sign(real(V(:,end)));



end

