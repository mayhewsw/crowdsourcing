function [ pHat, tHat ] = discretized_bp( A )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% dummy variable.
pHat = 0;

% find a list of all edges.
[row,col,v] = find(A~=0);
edges = [row col];

% note: first col of edges is tasks, second col is workers. Each tuple
% (i,j)
% represents an edge from task i to worker j.

alpha = 6;
beta = 2;

% randomly initialize x
xmap = containers.Map();

for r=1:size(edges,1)
    edge = edges(r,:);
    a = edge(1);
    b = edge(2);
    
    xmap(frmt(a,b,+1)) = rand;
    xmap(frmt(a,b,-1)) = rand;
    
end


% because why not.
T = 10;

% how many chunks to split p into?
chunks = 5;
pVals = 0.01:(1/chunks):1;
pVals(end) = 0.99; % just keeps everything in the range of 0,1...

ymap = containers.Map();

for iter=1:T
    fprintf('iteration %d\n', iter);
    % update y
    for p=pVals
        for r=1:size(edges,1);
            % reverse edge because it is y...
            edge = edges(r,:);
            i = edge(1);
            a = edge(2);
            
            % get neighbors of a (column in A)
            [nbrs,~,~] = find(A(:,a)~=0);
            
            % get neighbors of a, excluding i
            prodTerm = 1;
            for j=nbrs' % matlab only iterates through row vectors.
                if j == i
                    continue
                end
                
                prodTerm = prodTerm * ((1 + (2*p-1)*A(j, a))*xmap(frmt(j,a,+1)) + ...
                    (1 - (2*p-1)*A(j, a))*xmap(frmt(j,a,-1)));
                
                if(prodTerm==0)
                    fprintf('is zero...');
                elseif(isnan(prodTerm))
                    fprintf('whoops... isnan\n')
                elseif(isinf(prodTerm))
                    fprintf('whoops... if inf\n')
                elseif(isnan(betapdf(p, alpha,beta)*prodTerm))
                    fprintf('oh whoops...\n');
                end
                
            end
            
            ymap(frmt(a,i,p)) = betapdf(p, alpha,beta) * prodTerm;
            
        end
    end
    
    % update x
    for tval=[-1 1]
        for r=1:size(edges,1);
            edge = edges(r,:);
            i = edge(1);
            a = edge(2);
            
            % get neighbors of i (row in A)
            [~,nbrs,~] = find(A(i,:)~=0);
            prodTerm = 1;
            for b=nbrs
                if b == a
                    continue;
                end
                
                innersum = 0;
                for pb=pVals
                    innersum = innersum + ymap(frmt(b,i,pb)) * (pb*(A(i,b)==tval) + (1-pb)*(A(i,b)~=tval));
                end
                prodTerm = prodTerm * innersum;
                
                if(isnan(prodTerm))
                    fprintf('whoops... isnan\n')
                elseif(isinf(prodTerm))
                    fprintf('whoops... if inf\n')
                end
                
            end
            
            xmap(frmt(i,a,tval)) = prodTerm;
        end
    end
end

% compute decision values.
[n,m] = size(A);
xdec = zeros(n,2);
for tval=[-1 1]
    for i=1:n
        
        % get neighbors of i (row in A)
        [~,nbrs,~] = find(A(i,:)~=0);
        prodTerm = 1;
        for b=nbrs
            innersum = 0;
            for pb=pVals
                innersum = innersum + ymap(frmt(b,i,pb)) * (pb*(A(i,b)==tval) + (1-pb)*(A(i,b)~=tval));
            end
            
            prodTerm = prodTerm * innersum;
            
        end
        
        if tval == -1
            ind=1;
        else
            ind=2;
        end
        
        xdec(i,ind) = prodTerm;
        
    end
end

% nice.
tHat = sign(xdec(:,2) - xdec(:,1));


end

function [s]=frmt(i,j,v)
s = sprintf('%d,%d,%f',i,j,v);
end

