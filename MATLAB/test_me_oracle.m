%%%%%%%%%%%%%%%%%%%%%%%
% Example of SEM-ME Recovery algorithm, given a noisy version of the mixing
% matrix.
%%%%%%%%%%%%%%%%%%%%%%%

p = 5;          % Number of underlying variables
p_z = 4;        % Number of Z variables (i.e., measured with error)
pr_edge = 0.4;  % Prob. of edge connection between each pair of variables
n = 10;         % Number of observations
delta = 0.01;   % Variance of the added Gaussian noise

% Data generating process
[A, ~, measure_idx, W, ~] = generate_me(p, p_z, pr_edge,n);
disp('True adjacency matrix:')
disp(A)

% Add Gaussian noise to the support of W
m = size(W, 2);
noise = delta * randn(p,m) .* (W ~= 0); 
% Add Gaussian noise to each entry of W with probability 0.2
noise = noise + delta * randn(p, m) .* randsrc(p, m, [0, 1; 0.8, 0.2]);

W = W + noise;
print_W = 0;
if print_W
    disp('Noisy mixing matrix:')
    disp(W)
end

% Recovery
[A_full, ~]=semme_recovery(W, 0.05, measure_idx);
disp('Recovered adjacency matrix:')
disp(A_full)