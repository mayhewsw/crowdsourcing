% crowdsourcing.

%% Parameters Section

% number of instances to average over
instances = 10;
n=30;
m=50;
a = 0.3;
b = 0.95;

% max size l to run to
lmax = 20;

% used to store 2 runs: simplified BP, and Majority Voting
runs = zeros(2, lmax-1);

%% Algorithm Section
for l=2:lmax
    
    avgerror = 0;
    avgmverror = 0;
    
    % each iteration of this loop is a completely new setup of the problem.
    for dontcare = 1:instances
        
        % task labels
        t = sign( rand(n,1)-0.5 );
        
        % worker reliabilities
        p = a+(b-a)*rand(m,1);
        
        [A,E] = generate_graph(p,t,l);
       
        [~, tHat, T_bp] = simplified_bp(A);
        
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
title(sprintf('Average Error over %d instances, each with %d iterations. m=%d, n=%d', instances, T_bp, m, n))
xlabel('l');
ylabel('P(Error)');

