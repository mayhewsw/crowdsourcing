function [ pHat, tHat ] = iwmv( A )

% Taken from:
% http://arxiv.org/pdf/1411.4086v1.pdf
[n,m] = size(A);

% initialization
v = ones(m,1);
T = A~=0;
L=2; % number of labels?

y = zeros(n,1);
w = zeros(m,1);

Z = A';

for iter=1:10
   
    % the y loop
    for j=1:n
        
        neg = 0;
        pos = 0;
        for i=1:m
            neg = neg + v(i) * (Z(i,j) == -1);
            pos = pos + v(i) * (Z(i,j) == 1);
        end
                
        if pos > neg
            y(j) = 1;
        else
            y(j) = -1;
        end
    end
    
    % the w loop
    for i=1:m
        num = 0;
        for j=1:n
            num = num + Z(i,j) == y(j);
        end
        denom = sum(T(i,:));
        w(i) = num / denom;
    end
    w(isnan(w)) = 0 ;
    
    % the v loop
    v = L*w -1;
    
end

pHat = 1;
tHat = -y;

end

