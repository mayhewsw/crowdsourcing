% crowdsourcing.

%% Parameters Section

% number of instances to average over
instances = 100;
n=250;
m=250;
a = 0.3;
b = 0.95;

% max size l to run to
lmax = 15;

% used to store 2 runs: simplified BP, and Majority Voting
runs = zeros(3, lmax-1);

%% Algorithm Section
for l=2:lmax
    
    avgerror_bp = 0;
    avgerror_mv = 0;
    avgerror_pi = 0;
    
    % each iteration of this loop is a completely new setup of the problem.
    for dontcare = 1:instances
        
        % task labels
        t = sign( rand(n,1)-0.5 );
        
        % worker reliabilities
        p = a+(b-a)*rand(m,1);
        
        [A,E] = generate_graph(p,t,l);
       
        [~, tHat_bp, T_bp] = simplified_bp(A);
        T_bp = 1;
        error_bp = sum(tHat_bp ~= t) / n;
        avgerror_bp = avgerror_bp + error_bp;
        
        error_mv = sum(sign(sum(A, 2)) ~= t) / n;
        avgerror_mv = avgerror_mv + error_mv;
        
        %[~, tHat_pi] = power_iteration(A);
        tHat_pi = ones(m,1);
        error_pi = sum(tHat_pi ~= t) / n;
        avgerror_pi = avgerror_pi + error_pi;
   
    end
    
    runs(1, l-1) = avgerror_bp / instances;
    runs(2, l-1) = avgerror_mv / instances;
    runs(3, l-1) = avgerror_pi / instances;
    
    
    fprintf('l=%d, Avg BP error: %f, Avg mv error: %f\n', l,avgerror_bp / instances, avgerror_mv / instances);
    
end

semilogy(2:lmax, runs(1,:), '-or', 2:lmax, runs(2, :), '-db', 2:lmax, runs(3,:), '-dk');
legend('Simplified BP', 'Majority Voting', 'Power Iteration');
title(sprintf('Average Error over %d instances, each with %d iterations. m=%d, n=%d', instances, T_bp, m, n))
xlabel('l');
ylabel('P(Error)');

