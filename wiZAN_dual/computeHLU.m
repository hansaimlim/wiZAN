function [Ru, Ru_max] = computeHLU(real, predict)
% real is {0,1} vector
% predict is [0,1] vector


[t,ind] = sort(predict,'descend');
rsize = size(real,1);
Ru = 0;
for i = 1:rsize
    if real(ind(i)) == 1
        Ru = Ru + 1 / (2^((i-1)/4));
    end
end
[I] = find(real);
isize = size(I,1);
if isize == 0
    Ru_max = 0;
else
    Ru_max = 0;
    for i = 1:isize
        Ru_max = Ru_max + 1 / (2^((i-1)/4));
    end
end
end

