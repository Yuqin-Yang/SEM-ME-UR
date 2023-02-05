function [A_full, uniq]=semme_recovery(W, alpha, index_measurement)
%%%%%%%%%%%%%%%%%%%%%%%
% SEM-ME recovery algorithm based on Direct Ordered Grouping (DOG).
%
% Input:
%   W:                 Recovered mixing matrix.
%   alpha:             Threshold for pruning.
%   index_measurement: Index of the Z variables
%
% Output:
%   A_full: Recovered adjacency matrix (among observed variables).
%   uniq:   Indicating whether the recovery is unique.
%%%%%%%%%%%%%%%%%%%%%%%

p = size(W,1);
uniq = 1;
W = prune_W_me(W, index_measurement);
[row_order, col_order] = find_aog_me(W, 0.05, 0);

W = W(:, col_order);
product = cartesianProduct(row_order);
p_s = size(W,2);
edge_est = p * p_s; % suff. large

for ii=1:size(product,1)
    nuleaf = product(ii,:);
    uleaf = setdiff(1:p, nuleaf);

    W = W ./ W(sub2ind(size(W), nuleaf, 1:p_s));
    A_new = W(nuleaf, :);
    A_new = inv(A_new);
    B_new = W(uleaf,:) * A_new;
    edge_new = sum(sum(abs(A_new) >= alpha)) + sum(sum(abs(B_new) >= alpha));
    if edge_new <= edge_est
        A = eye(p_s)- A_new;
        A(abs(A) < alpha) = 0;
        B = B_new;
        B(abs(B) < alpha) = 0;
        A_full = zeros(p,p);
        A_full(nuleaf, nuleaf) = A;
        A_full(uleaf, nuleaf) = B;
        edge_est = edge_new;
        if (edge_new == edge_est) %&& min()
            uniq = 0;
        end
    end
end
end

function W = prune_W_me(W, index_measure)
%%%%%%%%%%%%%%%%%%%%%%%
% Extract the mixing matrix corresponding to the underlying model (i.e.,
% W^ME) from W, and prune the mixing matrix W by iteratively removing
% the entries with the smallest absolute value, until a valid AOG can be 
% recovered.
%%%%%%%%%%%%%%%%%%%%%%%

p = size(W, 1);
m = size(W, 2);

%%%% Remove the columns that corresponds to the measurement error
W = W ./ max(vecnorm(W), 1e-6);
for i=1:m
     [~,index]=max(abs(W(:,i)));
     if W(index,i)<0
         W(:,i)= -1*W(:,i);
     end
end
AA = eye(p);
AA = AA(:,index_measure);
Cost = zeros(length(index_measure),m);
    for i=1:m
        Cost(:,i) = sum(abs(AA - W(:,i)));
    end

col = matchpairs(Cost,1000);
W(:,col(:,2)) = [];

%%%% Pruning
W_copy = W;

[~, idx] = sort(abs(W(:)));
col_order = [];
i = max(sum(W(:) == 0), 1);
while length(col_order) < size(W, 2)
    W(idx(i)) = 0;
    [~, col_order] = find_aog_me(W, 0.05,0);
    i = i + 1;
    if i == length(W(:))  % another run with normalization
        W = W_copy;
        col = matchpairs(-W,1000);
        W = W ./ W(sub2ind(size(W),col(:,1),col(:,2)))';
        [~, idx] = sort(abs(W(:)));
        col_order = [];
        j = max(sum(W(:) == 0), 1);
        while length(col_order) < size(W, 2)
            W(idx(j)) = 0;
            [~, col_order] = find_aog_me(W, 0.05,0);
            j = j + 1;
        end
    end
end
end

function [row_order, col_order] = find_aog_me(W, alpha, leaf_node)
%%%%%%%%%%%%%%%%%%%%%%%
% Find the Ancestral Ordered Grouping (AOG) of the model (cf. Alg. 1).
%
% Input:
%   W:         Recovered mixing matrix.
%   alpha:     Threshold for pruning.
%   leaf_node: Indicating whether to return groups of only uleaf node.
%              If True then all groups are returned. Otherwise only groups
%              containing nuleaf nodes are returned.
%
% Output:
%   row_order: AOG of the model.
%   col_order: A Causal order among the observed variables.
%%%%%%%%%%%%%%%%%%%%%%%

W(abs(W) < alpha) = 0;
n_nz = sum(W ~= 0, 2);
p = size(W,1);
p_s = size(W,2); % number of sources, p_s < p

col_index = 1:p_s;
row_index = 1:p;
row_order = cell(p_s,1);
col_order = [];
leaf = [];

for jj=1:p_s
    w_temp = W(row_index, col_index);
    temp_index = find(sum(w_temp ~= 0, 2)==1);
    [n0, idx] = min(n_nz(temp_index));
    row = temp_index(idx);
    col = find(w_temp(row,:));
    temp_index = temp_index(w_temp(temp_index, col) ~= 0);
    non_leaf = [];
    
    for ii = 1:length(temp_index)
       temp_row = temp_index(ii);
       if n_nz(temp_row) == n0
           non_leaf = [non_leaf; row_index(temp_row)];
       else
           leaf = [leaf; row_index(temp_row)];
       end
    end
    row_order(jj) = {non_leaf};
    col_order = [col_order; col_index(col)];
    
    row_index(temp_index) = [];
    col_index(col) = [];
    n_nz(temp_index) = [];
    
end
if leaf_node && any(leaf) 
    row_order(p_s+1) = {leaf};
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