% this file creates UAI format configuration file for Dawid-Skene model 
% no of workers 
m = 5; 
% no of items 
n = 10; 
discretizationSize = 3; 

alpha = 6;
beta = 2;
% task labels
t = sign( rand(n,1)-0.75 );
l = 1; 
% worker reliabilities
p = ones(m, 1); % 0.1+0.9*betarnd(alpha,beta,m,1);
[A,E] = generate_graph(p,t,l);
E = min(E,1); % just to make sure the elements of E are not more than one 

% write into file 
id = fopen('output.uai', 'w+');
% header
fprintf(id, 'MARKOV\n');
% number of variables 
fprintf(id, '%d\n', n + m);
% number of values each variable can take on 
for i = 1:m
    fprintf(id, '%d ', discretizationSize);
end 
for i = 1:n
    fprintf(id, '%d ', 2);
end 

% number of factors 
fprintf(id, '\n%d\n', sum(sum(E)));

% specification of each factor 
for i = 1:n 
    for j = 1:m 
        if A(i, j) ~= 0 
            fprintf(id, '2 %d %d \n', j-1, m + i-1);  % -1 for starting from zero 
        end 
    end
end

% probability tables for each factor 
for i = 1:n 
    for j = 1:m 
        if A(i, j) == 0 
            continue; 
        end 
        fprintf(id, '%d\n', 2*discretizationSize);
        for workerValue = 1:discretizationSize
            for label = [-1, 1]
                if A(i, j) == label 
                    fprintf(id, '%f ', 1.0 * ( workerValue - .5 ) / discretizationSize);
                else
                    fprintf(id, '%f ', 1 - 1.0 * ( workerValue - .5 ) / discretizationSize);                    
                end
            end
            fprintf(id, '\n');
        end 
    end 
end
fclose(id)
