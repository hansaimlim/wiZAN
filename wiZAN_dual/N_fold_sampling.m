function N_fold_sampling(R_whole, Sample_mat, N, outdir)
%returns training and test sets 10-folded from Sample_mat (could be same as R_whole if whole sample test)
%Sample_mat must be a part of the R_whole matrix (could be equal)
%Each training set = R_whole - test_set
%N=10;	%N fold
numedge_whole=sum(R_whole(:)>0);
numedge_sample=sum(Sample_mat(:)>0);
siz=size(R_whole);
mkdir(outdir);

b=floor(numedge_sample/N);	%base num edges in testfile
r=rem(numedge_sample, N);	%remainder
%testarray=createArrays(N, siz);	%cell array to contain test sets
testcumulative=zeros(size(Sample_mat));	%to check sum(all_tests) == Sample_mat
Sample_up=Sample_mat; %Sample_up will be updated each iteration
for n=1:N
    num=b;	%num edge in test
    if r>0
        num = num + 1;
        r = r - 1;
    end
    idx = find(Sample_up);	%must sample from updated R
    idxpick = datasample(idx, num, 'Replace', false);	%pick num random indexes from idx
    [i j] = ind2sub(siz, idxpick);
    test = sparse(i(:), j(:), 1, siz(1), siz(2));
    testcumulative = testcumulative + test;
    train = R_whole - test;	%train + test = R
    Sample_up = Sample_up - test;	%update R
    [z x c] = find(test);
    [j k l] = find(train);
    testfile=[outdir 'test' num2str(n) '.csv'];
    trainfile=[outdir 'train' num2str(n) '.csv'];
    csvwrite(testfile, [z, x]);
    csvwrite(trainfile, [j, k]);
end

%to compare the sum of all test sets with Sample_mat
mustequal = isequal(testcumulative, Sample_mat);
if mustequal
 disp('N-fold sampling well done!')
else
 msg = 'Sum of all test sets NOT equal to the given sample matrix. Please check the code.';
 error(msg);
end

function result = createArrays(nArrays, arraySize)
    result = cell(1, nArrays);
    for i = 1 : nArrays
        result{i} = zeros(arraySize);
    end
end

end
