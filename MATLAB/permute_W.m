function [B, B_learn, error, struc] = permute_W(B, B_learn)
%%%%%%%%%%%%%%%%%%%%%%%
% Repermute and rescale B_learn such that its columns has the same order 
% and scale as B.
%
% Input:
%   B:          True matrix. 
%   B_learn:    Recovered matrix.
%
% Output:
%   B:          True matrix (after rescaling). 
%   B_learn:    Recovered matrix (after repermutation and rescaling).
%   error:      Frobenious norm between B and B_learn.
%   struc:      Indicating whether B and B_learn has the same support.
%%%%%%%%%%%%%%%%%%%%%%%

m = size(B,2);
B = B ./ vecnorm(B);
for i=1:m
     [~,index]=max(abs(B(:,i)));
     if B(index,i)<0
         B(:,i)= -1*B(:,i);
     end
end

B_learn = B_learn ./ max(vecnorm(B_learn), 1e-6);
for i=1:m
     [~,index]=max(abs(B_learn(:,i)));
     if B_learn(index,i)<0
         B_learn(:,i)= -1*B_learn(:,i);
     end
end

Cost = zeros(m,m);
    for i=1:m
        Cost(i,:) = sum(abs(B - B_learn(:,i)));
    end
    
col = matchpairs(Cost,1000);
B_learn = B_learn(:,col(:,1));
error = norm(B - B_learn);
struc = sum(sum(xor(B~=0, B_learn~=0)));



