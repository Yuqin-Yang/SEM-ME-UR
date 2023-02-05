%%%%%%%%%%%%%%%%%%%%%%%
% Example of SEM-ME Recovery algorithm, given observed data generated from 
% a model with non-Gaussian noise.
%%%%%%%%%%%%%%%%%%%%%%%

p = 4;                          % Number of underlying variables
p_znl = 1;                      % Number of Z^NL variables
pr_edge = min(0.5, 2.5/(p-1));  % Prob. of edge connection
n = 500 * (p + p_znl);          % Number of samples

% Data generating process
[A, ~, measure_idx, W, X] = generate_me_znl(p, p_znl, pr_edge,n);
disp('True adjacency matrix:')
disp(A)
    
% Recover W using Overcomplete ICA
m = size(W, 2);
W_learn = [];
while ~any(any(W_learn)) % W_learn do not have empty columns
    W_learn = support(X', m, 0.8*n, 50, 0.2);
end
[W, W_learn, ~, ~] = permute_W(W, W_learn);
disp('True mixing matrix:')
disp(W)
disp('Recovered mixing matrix (repermuted according to true W):')
disp(W_learn)

% Recovery
[A_full, ~]=semme_recovery(W_learn, 0.05, measure_idx);
disp('Recovered adjacency matrix:')
disp(A_full)