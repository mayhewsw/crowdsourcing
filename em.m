function [ pHat, tHat ] = em( A )

% Taken from:
% http://web.engr.illinois.edu/~swoh/courses/IE598/handout/crowd.pdf
[n,m] = size(A);

% num pos responses / nbrs of i
q = zeros(n,2);
q(:,1) = sum(A==1, 2) ./ sum(A~=0, 2);
q(isnan(q)) = 0 ;

q(:,2) = 1 - q(:,1);

% initialize pHat
pHat = ones(m,1);

% Horrible matlab, I know, but done is better than beautiful.
for j=1:m
    s = 0;
    tot = 0;
    for i=1:n
        if A(i,j) ~= 0
            if A(i,j) == 1
                s = s + q(i,1);
            elseif A(i,j) == -1
                s = s + q(i,2);
            end
            tot = tot + 1;
        end
    end
    pHat(j) = s / tot;
end

pHat(isnan(pHat)) = 0 ;


for iter=1:10
    
    % E step
    for i=1:n
        numposprod = 1;
        numnegprod = 1;
        
        for j=1:m
            % only nbrs of i
            if A(i,j) ~= 0
                numposprod = numposprod * (pHat(j) * (A(i,j) == 1) + (1-pHat(j)) * (A(i,j) ~= 1));
                numnegprod = numnegprod * (pHat(j) * (A(i,j) == -1) + (1-pHat(j)) * (A(i,j) ~= -1));
            end
        end
        
        denom = numposprod + numnegprod;
        q(i, 1) = numposprod / denom;
        q(i, 2) = numnegprod / denom;
    end
    
    q(isnan(q)) = 0 ;
    
    % M step
    % Horrible matlab, I know, but done is better than beautiful.
    for j=1:m
        s = 0;
        tot = 0;
        for i=1:n
            if A(i,j) ~= 0
                if A(i,j) == 1
                    s = s + q(i,1);
                elseif A(i,j) == -1
                    s = s + q(i,2);
                end
                tot = tot + 1;
            end
        end
        pHat(j) = s / tot;
    end
    
    pHat(isnan(pHat)) = 0 ;
    
end

% get the max indices for each row
[~, tHat] = max(q, [], 2);

% convert these (1s and 2s) to -1s and 1s respectively;
% (note that this is backwards: the first column is P(t=+1) and the 
% second column is P(t=-1))
tHat = (tHat - 1)*2 - 1;

% since it is backwards, reverse it.
tHat = -tHat;
end

