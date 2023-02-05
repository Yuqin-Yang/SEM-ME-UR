function [A, B, uniq]=semur_recovery(W, alpha)
%%%%%%%%%%%%%%%%%%%%%%%
% SEM-UR recovery algorithm based on Direct Ordered Grouping (DOG).
%
% Input:
%   W:     Recovered mixing matrix.
%   alpha: Threshold for pruning.
%
% Output:
%   A:     Recovered adjacency matrix (among observed variables).
%   B:     Recovered latent connection matrix.
%   uniq:  Indicating whether the recovery is unique.
%%%%%%%%%%%%%%%%%%%%%%%

p = size(W,1);
p_s = size(W,2);
uniq = 1;
W = prune_aog(W);
[row_order, col_order] = find_aog(W, 0.05, 0);

W = W(row_order,:);
product = cartesianProduct(col_order);

edge_est = p * p_s; % suff. large
A = zeros(p,p);
B = zeros(p, p_s - p);

for ii=1:size(product,1)
    observed = product(ii,:);
    latent = setdiff(1:p_s, observed);

    A_new = W(:, observed); % (I - A)^-1
    A_new = A_new ./ diag(A_new)';
    A_new = inv(A_new);
    B_new = A_new * W(:, latent);
    edge_new = sum(sum(abs(A_new) >= alpha)) + sum(sum(abs(B_new) >= alpha));
    if edge_new <= edge_est
        A = eye(p)- A_new;
        A(abs(A) < alpha) = 0;
        B = B_new;
        B(abs(B) < alpha) = 0;
        edge_est = edge_new;
        if (edge_new == edge_est)
            uniq = 0;
        end
    end
end
rev(row_order) = 1:p;
A = A(rev, rev);
B = B(rev,:);
end

function W = prune_aog(W)
%%%%%%%%%%%%%%%%%%%%%%%
% Prune the mixing matrix W by iteratively removing the entries with the 
% smallest absolute value, until a valid AOG can be recovered.
%%%%%%%%%%%%%%%%%%%%%%%

[~, idx] = sort(abs(W(:)));
row_order = [];
i = max(sum(W(:) == 0), 1);
while length(row_order) < size(W, 1)
    W(idx(i)) = 0;
    [row_order, ~] = find_aog(W, 0.05,0);
    i = i + 1;
end
end

function [row_order, col_order] = find_aog(W, alpha, leaf_node)
%%%%%%%%%%%%%%%%%%%%%%%
% Find the Ancestral Ordered Grouping (AOG) of the model (cf. Alg. 1).
%
% Input:
%   W:         Recovered mixing matrix.
%   alpha:     Threshold for pruning.
%   leaf_node: Indicating whether to return groups of only latent variable.
%              If True then all groups are returned. Otherwise only groups
%              containing observed variables are returned.
%
% Output:
%   row_order: A Causal order among the observed variables.
%   col_order: AOG of the model.
%%%%%%%%%%%%%%%%%%%%%%%

W(abs(W) < alpha) = 0;
n_nz = sum(W ~= 0, 1);
p = size(W,1);
p_s = size(W,2); % number of sources

column_index = 1:p_s;
row_index = 1:p;
col_order = cell(p,1);
row_order = [];
leaf = [];

for jj=1:p
    w_temp = W(row_index, column_index);
    temp_index = find(sum(w_temp ~= 0, 1)==1);
    [n0, idx] = min(n_nz(temp_index));
    col = temp_index(idx);
    row = find(w_temp(:,col));
    temp_index = temp_index(w_temp(row, temp_index) ~= 0);
    non_leaf = [];
    
    for ii = 1:length(temp_index)
       temp_col = temp_index(ii);
       if n_nz(temp_col) == n0
           non_leaf = [non_leaf; column_index(temp_col)];
       else
           leaf = [leaf; column_index(temp_col)];
       end
    end
    col_order(jj) = {non_leaf};
    row_order = [row_order; row_index(row)];
    
    column_index(temp_index) = [];
    row_index(row) = [];
    n_nz(temp_index) = [];
    
end
if leaf_node && any(leaf) 
    col_order(p+1) = {leaf};
end
end

function result = cartesianProduct(sets)
%%%%%%%%%%%%%%%%%%%%%%%
% Return the Cartesian Product of the elements in the set.
%%%%%%%%%%%%%%%%%%%%%%%

c = cell(1, numel(sets));
[c{:}] = ndgrid( sets{:} );
result = cell2mat( cellfun(@(v)v(:), c, 'UniformOutput',false) );
end