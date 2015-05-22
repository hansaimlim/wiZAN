function [PR_1, PR_2] = computePR(real, predict)
% real is {0,1} vector
% predict is [0,1] vector


[t,ind] = sort(predict,'descend');
rsize = size(real,1);
PR_1 = 0;
for i = 1:rsize
    if real(ind(i)) == 1
        PR_1 = PR_1 + (i-1)/rsize;
    end
end
[I] = find(real);
PR_2 = size(I,1); 
end

