% crowdsourcing.
close all
clear

%% Parameters Section

% number of instances to average over
instances = 3;
n=100;
m=100;
alpha = 6;
beta = 2;

% max size l to run to
lmax = 15;

algs_to_run = {'mv', 'discretized_bp'};
numalgs = numel(algs_to_run);

% used to store all runs.
runs = zeros(numalgs, lmax-1);

%% Algorithm Section
for l=15:lmax
    fprintf('l=%d\n', l);
    
    allruns = zeros(numalgs, instances);
    
    % each iteration of this loop is a completely new setup of the problem.
    for inst = 1:instances
        
        % task labels
        t = sign( rand(n,1)-0.75 );
        
        % worker reliabilities
        p = 0.1+0.9*betarnd(alpha,beta,m,1);
        
        [A,E] = generate_graph(p,t,l);
        
        for i=1:numalgs
            alg = algs_to_run(i);
            [~, tHat] = eval(sprintf('%s(A)', alg{1}));
            error = sum(tHat ~= t) / n;
            
            allruns(i,inst) = error;
        end
    end
    
    means = mean(allruns, 2);
    stds = std(allruns, 0, 2) / sqrt(instances);
    
    scores = means + stds
    
    runs(:, l-1) = scores;
    
    
    
end

cc=hsv(numalgs*3);

for i=1:numalgs
    semilogy(2:lmax, runs(i,:), 'color',cc(i*3,:));
    hold on
end


legend(algs_to_run, 'Interpreter', 'none');
title(sprintf('Average Error over %d instances. m=%d, n=%d', instances, m, n))
xlabel('l');
ylabel('P(Error)');

annotation('textbox', [0.14,0.14,0.14,0.05],...
    'String', sprintf('best accuracy at l=15: %f', min(runs(:,end))));


