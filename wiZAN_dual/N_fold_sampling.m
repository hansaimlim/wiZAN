function N_fold_sampling(N, R_csv, nrow, ncol, outdir)
%returns training and test sets from input R matrix
%each training+test equals to the R (e.g. train1 + test1 = train2 + test2 = R)
%N=10;	%N fold
%rline='some_csv_file.csv';
%nrow=198712;
%ncol=3549;
rline=csvread(R_csv);
R=sparse(rline(:,1), rline(:,2), 1, nrow, ncol);
numedge=sum(R(:)>0);
siz=size(R);

Rupdate=R;	%matrix R updated by iteration
b=floor(numedge/N);	%base num edges in testfile
r=rem(numedge, N);	%remainder
%testarray=createArrays(N, siz);	%cell array to contain test sets
testcumulative=zeros(nrow,ncol);	%to check sum(all_tests) == R
for n=1:N
 num=b;	%num edge in test
 if r>0
  num = num + 1;
  r = r - 1;
 end
idx = find(Rupdate);	%must sample from updated R
idxpick = datasample(idx, num, 'Replace', false);	%pick num random indexes from idx
[i j] = ind2sub(siz, idxpick);
test = sparse(i(:), j(:), 1, siz(1), siz(2));
testcumulative = testcumulative + test;
train = R - test;	%train + test = R
Rupdate = Rupdate - test;	%update R
[z x c] = find(test);
[j k l] = find(train);
testfile=[outdir 'test' num2str(n) '.csv'];
trainfile=[outdir 'train' num2str(n) '.csv'];
csvwrite(testfile, [z, x]);
csvwrite(trainfile, [j, k]);
end

%to compare the sum of all test sets with R
mustequal = isequal(testcumulative, R);
if mustequal
 mustequal
else
 msg = 'Sum of all test sets NOT equal to R. Please check the code.';
 error(msg);
end

function result = createArrays(nArrays, arraySize)
    result = cell(1, nArrays);
    for i = 1 : nArrays
        result{i} = zeros(arraySize);
    end
end

end
