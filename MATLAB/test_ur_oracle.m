%%%%%%%%%%%%%%%%%%%%%%%
% Example of SEM-UR Recovery algorithm, given a noisy version of the mixing
% matrix.
%%%%%%%%%%%%%%%%%%%%%%%

q = 4;        % Number of observed variables
l = 1;        % Number of latent variables
pr_obs = 0.4; % Prob. of edge connection between obs. variables
pr_lat = 0.4; % Prob. of edge connection between lat. & obs. variable
n = 10;       % Number of observations
delta = 0.01;   % Variance of the added Gaussian noise.

% Data generating process
[A,B,W,~] = generate_ur(q, l, pr_obs, pr_lat,n);
disp('True adjacency matrix (A & B):')
disp(A)
disp(B)

% Add Gaussian noise to the support of W
m = size(W, 2);
noise = delta * randn(q, m) .* (W ~= 0); 
% Add Gaussian noise to each entry of W with probability 0.2
noise = noise + delta * randn(q, m) .* randsrc(q, m, [0, 1; 0.8, 0.2]);

W = W + noise;
print_W = 0;
if print_W
    disp('Noisy mixing matrix:')
    disp(W)
end

% Recovery
[A_rec, B_rec, ~]=semur_recovery(W, 0.05);
disp('Recovered adjacency matrix (A & B):')
disp(A_rec)
disp(B_rec)