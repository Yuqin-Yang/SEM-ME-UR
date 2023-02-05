function [A, order, measure_idx, W, X] = generate_me_znl(p, p_znl, pr_edge,n) 
%%%%%%%%%%%%%%%%%%%%%%%
% Data generating process of a linear SEM-ME, given the number of 
%   unobserved non-leaf variables (Z^NL).
% Designed for synthetic data simulation; the mixing matrix is 
%   of dimension p * (p + p_znl).
%
% Input:
%   p:       Number of underlying variables.
%   p_znl:   Number of unobserved non-leaf variables (Z^NL). Notice that
%            nu-leaf nodes include Z^NL and observed variables (Y). 
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
% Ensures that there are at least p_znl non-leaves in the model.
A = 0;
while sum(any(A)) < p_znl
    A = randsrc(p,p,[0,1;1-pr_edge,pr_edge]);
    A = tril(A,-1);
end

A = A .* (0.5 + 0.5 * rand(p,p));
W = inv(eye(p)-A);
perm_row = randperm(p);
A = A(perm_row, perm_row);
W = W(perm_row, :);
order(perm_row) = 1:p;

%%%%  Add measurement error to the mixing matrix  %%%%
leaf = 1:p;
leaf = leaf(all(~A));
nl = setdiff(1:p, leaf);
measure_ul = randsample(nl, p_znl, false);
measure_idx = measure_ul;
% Each leaf node is observed with probability 0.5
for i=leaf
    if rand() < 0.5
        measure_idx = [measure_idx, i];
    end
end

mat = eye(p);
W = [W, mat(:,measure_ul)];
m = size(W, 2);
perm_col = randperm(m);
W = W(:, perm_col);

%%%%  Non-gaussian noise  %%%%
% Noises are randomly drawn from Uniform[-0.5, 0.5]
X = W * (rand(m, n) - 0.5);
end