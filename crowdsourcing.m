% crowdsourcing.

% number of instances to average over
instances = 10;
n=100;
m=100;
a = 0.3;
b = 0.95;

% max number of iterations
T = 10;

% max size l to run to
lmax = 20;

runs = zeros(2, lmax-1);

for l=2:lmax
    
    avgerror = 0;
    avgmverror = 0;
    
    for dontcare = 1:instances
        
        % task labels
        t = sign( rand(n,1)-0.5 );
        
        % worker reliabilities
        p = a+(b-a)*rand(m,1);
        
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
     
        % x,y is initialized for every edge i,j
        x = sparse(zeros(n,m));
        y = sparse(normrnd(1,1,m,n));
        
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
        
        errorrate = sum(tHat ~= t) / n;
        avgerror = avgerror + errorrate;
        
        mverror = sum(sign(sum(A, 2)) ~= t) / n;
        avgmverror = avgmverror + mverror;
   
    end
    
    runs(1, l-1) = avgerror / instances;
    runs(2, l-1) = avgmverror / instances;
    
    fprintf('l=%d, Avg BP error: %f, Avg mv error: %f\n', l,avgerror / instances, avgmverror / instances);
    
end

plot(2:lmax, runs(1,:), '-or', 2:lmax, runs(2, :), '-db');
legend('Simplified BP', 'Majority Voting');
title(sprintf('Average Error over %d instances, each with %d iterations. m=n=%d', instances, T, n))
xlabel('l');
ylabel('P(Error)');

