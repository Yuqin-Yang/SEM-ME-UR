function [A, order, measure_idx, W, X] = generate_me(p, p_z, pr_edge,n) 
%%%%%%%%%%%%%%%%%%%%%%%
% Data generating process of a linear SEM-ME, given the number of 
% underlying variables.
%
% Input:
%   p:       Number of underlying variables.
%   p_z:     Number of Z variables (i.e., measured with error).
%   pr_edge: Prob. of edge connection between each pair of variables.
%   n:       Number of observations (samples/measurements) for each 
%            observed variable.
%
% Output:
%   A:           True adjacency matrix (among underlying variables).
%   order:       Causal order among variables.
%   measure_idx: Index of the variables that are measured with error.
%   W:           True mixing matrix (not W^ME).
%   X:           Observed data.
%%%%%%%%%%%%%%%%%%%%%%%

%%%%  Generate matrix A  %%%%
% Ensures that the adj. matrix A is non-zero
A = 0;
while ~any(any(A))
    A = randsrc(p,p,[0,1;1-pr_edge,pr_edge]);
    A = tril(A,-1);
end

A = A .* (0.5 + 0.5 * rand(p,p));
W = (eye(p)-A) \ eye(p);
perm_row = randperm(p);
A = A(perm_row, perm_row);
W = W(perm_row, :);
order(perm_row) = 1:p;

%%%%  Add measurement error to the mixing matrix  %%%%
measure_idx = randsample(p, p_z, false);
leaf = 1:p;
leaf = leaf(all(~A));
nu_leaf = setdiff(measure_idx, leaf);
mat = eye(p);
W = [W, mat(:,nu_leaf)];
m = size(W, 2);
perm_col = randperm(m);
W = W(:, perm_col);


%%%%  Non-gaussian noise  %%%%
% Noises are randomly drawn from Uniform[-0.5, 0.5]

X = W * (rand(m, n) - 0.5);
end