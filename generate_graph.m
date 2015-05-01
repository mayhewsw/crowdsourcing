function [ A,E ] = generate_graph( p,t,l )
% Generate a graph
%   This implementation is erdos-renyi

m = numel(p);
n = numel(t);

% generate random graph
E = ceil( rand(n,m)-1+(l/m) );

A = rand(n,m) .* E;

pmat = repmat(p', n, 1);

% these indices will get t(i)
posT = A<pmat & A~=0;

% these indices will get -t(i)
negT = A>=pmat & A~=0;

tmat = repmat(t, 1, m);

A(posT) = tmat(posT);
A(negT) = -tmat(negT);


end

