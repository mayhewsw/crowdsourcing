function [ A,E ] = generate_graph( p,t,l )
% Generate a graph
%   This implementation is erdos-renyi

m = numel(p);
n = numel(t);

% generate random graph
E = ceil( rand(n,m)-1+(l/m) );
A = zeros(n,m);

% probably a more efficient way to do this...
for i=1:n
    for j=1:m
        if E(i,j) == 1
            if rand < p(j)
                A(i,j) = t(i);
            else
                A(i,j) = -t(i);
            end
        end
    end
end


end

