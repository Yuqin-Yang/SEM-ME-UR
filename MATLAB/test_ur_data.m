%%%%%%%%%%%%%%%%%%%%%%%
% Example of SEM-UR Recovery algorithm, given observed data generated from 
% a model with non-Gaussian noise.
%%%%%%%%%%%%%%%%%%%%%%%

q = 4;        % Number of observed variables
l = 1;        % Number of latent variables
pr_edge = min(0.5, 2.5/(q-1));  % Prob. of edge connection
n = 500 * (q + l);          % Number of samples

% Data generating process
[A, B, W, X] = generate_ur(q, l, pr_edge, pr_edge,n);
disp('True adjacency matrix (A & B):')
disp(A)
disp(B)
    
% Recover W using Overcomplete ICA
m = size(W, 2);
W_learn = [];
while ~any(any(W_learn)) % W_learn do not have empty columns
    W_learn = support(X', m, 0.8*n, 50, 0.1);
end
[W, W_learn, ~, ~] = permute_W(W, W_learn);
disp('True mixing matrix:')
disp(W)
disp('Recovered mixing matrix (repermuted according to true W):')
disp(W_learn)

% Recovery
[A_rec, B_rec, ~]=semur_recovery(W_learn, 0.05);
disp('Recovered adjacency matrix (A & B):')
disp(A_rec)
disp(B_rec)