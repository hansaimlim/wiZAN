function AP = computeAP(real, predict)
% real is {0,1} vector
% predict is [0,1] vector
% topN is 10


[t,ind] = sort(predict,'descend');
rsize = size(real,1);
prec = 0;
prec_sum = 0;
for i = 1:rsize
    if real(ind(i)) == 1
        prec = prec + 1;
        prec_sum = prec_sum + prec / i;
    end
end
[I] = find(real);
isize = size(I,1);
if isize == 0
    AP = 0;
else
    AP = prec_sum / isize;
end
end

