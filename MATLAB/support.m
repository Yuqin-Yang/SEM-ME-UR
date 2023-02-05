function Bout = support(V, m, train_size, bootn, th)
%%%%%%%%%%%%%%%%%%%%%%%
% Recover the mixing matrix using Overcomplete ICA and bootstrapping.
%
% Input:
%   V:          Observed data ([# of samples] * [# of variables]). 
%   m:          Number of sources (assumed to be known).
%   train_size: Number of samples in each bootstrapping iteration.
%   bootn:      Number of bootstrapping iterations.
%   th:         Pruning threshold.
%
% Output:
%   Bout:       Recovered mixing matrix.
%%%%%%%%%%%%%%%%%%%%%%%

p = size(V, 2);
check = 0;
while ~check
    Bhbase = estimate_model(V, m, train_size, 1);
    check = (size(Bhbase, 2) == m);
end

Btot = zeros(p, m, bootn);
Btot(:,:,1) = Bhbase;
for k=2:bootn
    check = 0;
    while ~check
        Bhnorm = estimate_model(V, m, train_size, 1);
        check = (size(Bhnorm, 2) == m);
    end   
    Cost = zeros(m,m);
    for i=1:m
        Cost(i,:) = sum((Bhbase-Bhnorm(:,i)).^2);
%         Cost(i,:) = sum(abs(Bhbase-Bhnorm(:,i)));
    end
    
    col = matchpairs(Cost,1000);
    Btot(:,:,k) = Bhnorm(:,col(:,1));
end

delta = 0.95;
nu = bootn - 1;
Bout = zeros(p,m);
Mu = mean(Btot, 3);
Std = zeros(p,m);
for i=1:p
    for j=1:m
        temp=reshape(Btot(i,j,:),1,bootn);
        s = std(temp);
        Std(i,j) = th + tinv((1+delta)/2,nu)*s/sqrt(bootn);
        if abs(Mu(i,j)) > Std(i,j)
            Bout(i,j) = Mu(i,j);
        end
    end
end
Bout = Bout./ max(abs(Bout));
end

function Bhnorm = estimate_model(X, m, train_size, replacement)
n = size(X,1);

C = cov(X);
D = 0;

if replacement == 1
    train_index = randsrc(1,train_size, [1:n;ones(1,n)/n]);
else
    train_index = randsample(n,train_size);
end

%%%%%%%%%%%%%%%% Reconstruction ICA %%%%%%%%%%%%%
if det(C) ~= 0
    [V, D] = eig(C);
    X = X * V * D^(-0.5); % prewhitening
end

Mdl = rica(X(train_index,:),m,'NonGaussianityIndicator',-1*ones(m,1));
W = Mdl.TransformWeights;

if any(D)
    W = V * D^(-0.5) * W;
end
W = pinv(W');
Bhnorm = W./ max(abs(W));

 for i=1:size(Bhnorm,2)
     [~,index]=max(abs(Bhnorm(:,i)));
     if Bhnorm(index,i)<0
         Bhnorm(:,i)= -1*Bhnorm(:,i);
     end
 end
end



