function [A,B,W,X] = generate_ur(q, l, pr_obs, pr_lat, n) 
%%%%%%%%%%%%%%%%%%%%%%%
% Data generating process of a linear SEM-UR.
%
% Input:
%   q:      Number of observed variables.
%   l:      Number of latent variables.
%   pr_obs: Prob. of edge connection between each pair of obs. variables.
%   pr_lat: Prob. of edge connection between each pair of 
%           latent variable and obs. variable.
%   n:      Number of observations (samples) for each observed variable.
%
% Output:
%   A:      True adjacency matrix (among observed variables).
%   B:      True latent connection matrix.
%   W:      True mixing matrix.
%   X:      Observed data.
%%%%%%%%%%%%%%%%%%%%%%%

%%%%  Generate matrix A  %%%%
% Ensures that the adj. matrix A is non-zero
A = 0;
while ~any(any(A))
    A = randsrc(q,q,[0,1;1-pr_obs,pr_obs]);
    A = tril(A,-1);
end
% Edge coefficients are randomly selected between 0.5 and 1
A = A .* (0.5 + 0.5 * rand(q,q));

%%%%  Generate matrix B  %%%%
% Ensures that each column of B has at least 2 non-zero entries

select = 0;
B = [];
while sum(select) < l
    B = randsrc(q,l,[0,1;1-pr_lat,pr_lat]);
    select = sum(B~=0) > 1;
end

if ~isempty(B)
    B = B .* (0.5 + 0.5 * rand(q,l));
end

B_full = [B eye(q)];
m = l + q;

W = (eye(q) - A) \ B_full;

perm_row = randperm(q);
perm_col = randperm(m);

W = W(perm_row, perm_col);
A = A(perm_row, perm_row);
B = B(perm_row, :);

%%%%  Non-gaussian noise  %%%%
% Noises are randomly drawn from Uniform[-0.5, 0.5]

X = W * (rand(m, n) - 0.5);
end