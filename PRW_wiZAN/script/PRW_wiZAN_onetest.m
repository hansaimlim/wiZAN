function [MAP, MPR, HLU, AUC] = PRW_wiZAN_onetest(train_csv, test_csv, prw_csv, outfile)
%modified version of wiZAN_dual for PRW_wiZAN process
%takes 3rd input argument 'prw' which contains prw result for all unique chemicals against all proteins in training set
%for single test pair (train, test, prw)
%fixed parameter as of 5/27/2015
%default rank=300, p6=0.75, p7=0.1 modified on 5/19/2015
para = [0.1, 300, 100, 0.75, 0.1]; % para: alpha, rank, maxIte, gamma, lambda

%pathways for necessary functions
addpath /scratch/hansaim.lim/wiZAN/wiZAN_dual/;
%chem_chem_zinc and protein_protein_zinc_blast matrices from chem-chem and prot-prot files
load /scratch/hansaim.lim/wiZAN/ZINC_data/chem_chem/chem_chem_zinc;
load /scratch/hansaim.lim/wiZAN/ZINC_data/prot_prot/protein_protein_zinc_blast;
%get number of chemical and protein
temp_c=size(chem_chem_zinc);
temp_p=size(protein_protein_zinc_blast);
m = temp_c(1);
n = temp_p(1);
%convert csv to matrix
train_line = csvread(train_csv);
train = sparse(train_line(:,1), train_line(:,2), 1, m, n);	%12384 chemicals and 3500 proteins in ZINC
prw_line = csvread(prw_csv);
P=sparse(prw_line(:,1),prw_line(:,2),prw_line(:,3), m, n);   %PRW results for test chemicals. Others are zeros
test = csvread(test_csv);

%protein_protein_zinc_blast = ceil(protein_protein_zinc_blast);
chem_chem_zinc = chem_chem_zinc + chem_chem_zinc';
%protein_protein_zinc_blast = protein_protein_zinc_blast + protein_protein_zinc_blast';
Pu=P.*(~train); %imPutation matrix. P(i, j) for test chemicals (only pairs NOT known). If a pair is known, value is 0. All others are zeros.
W=train+Pu; %weight, 1 for train pairs

summ = sum(chem_chem_zinc,2); %sum by rows
Dm = spdiags(summ,0,m,m);
Lu = Dm - chem_chem_zinc;

sumn = sum(protein_protein_zinc_blast,2); %sum by rows
Dn = spdiags(sumn,0,n,n);
Lv = Dn - protein_protein_zinc_blast;

[U, V] = updateUV(train, Lu, Lv, para, W, P);

%get predicted scores and ranks based on updated U and V
test_result = TPRbyRowRank(get_test_result(test, U, V), 200);	%max cutoff rank 200

%[MPR_C, MPR_U, MPR_I] = get_diff_coldstart(train, test, U, V);
[MAP, MPR, HLU, AUC] = get_diff(test, U, V, para);
fprintf('MAP = %0.4f, MPR = %0.4f, HLU = %0.4f, AUC = %0.4f\n', MAP, MPR, HLU, AUC);
outfileId=fopen(outfile, 'w');
fprintf(outfileId, '%5d %12.8f\n', test_result');
fclose(outfileId);
fprintf('Result file saved: %s\n',outfile);
clear train;
clear test;
clear test_result;
clear P;
clear para;
end

function test_result = get_test_result(test, U, V)
%to get the actual scores (and rank)
%result matrix contains 5 columns: chemical index, protein index, predicted
%score, row rank, column rank
scores = U * V';
test_size = size(test);
test_row_size = test_size(1);
test_result = zeros(test_row_size, 5);
for test_r = 1:test_row_size
    chem_index = test(test_r, 1);
    prot_index = test(test_r, 2);
    pred_score = scores(chem_index, prot_index);
    
    temp_row = scores(chem_index, :);
    row_rank = sum(temp_row >= pred_score) + 1; %the dense rank of the prediction score in the row (chemicals)
    temp_col = scores(:, prot_index);
    col_rank = sum(temp_col >= pred_score) + 1; %the rank of the prediction score in the column (proteins)
    test_result(test_r, :) = [chem_index, prot_index, pred_score, row_rank, col_rank];
end
end

function TPRbyRow = TPRbyRowRank(testResult, maxRank)
%to get TPR by cutoff Row rank
%output contains vector of dimension (maxRank, 2)
%testResult input must contain row rank on the 4th column
%testResult input must contain ONLY test pairs, since the number of lines will be used as condition positive
[cp, col] = size(testResult); %cp: condition positive, col: number of columns
TPRbyRow = zeros(maxRank, 2);
row_rank = testResult(:, 4);	%Rank by Rows on 4th column
rank = 1;	%cutoff rank
while rank <= maxRank
	tpr = sum(row_rank<=rank)/cp;
	TPRbyRow(rank, :) = [rank, tpr];	%update result row
	rank = rank + 1;
end
end

function [U, V] = updateUV(R, Lu, Lv, para, W, P)
% para: lambda, r, T, rank, maxIte, ite_of_bisection method, topN
[m, n] = size(R);
alpha = para(1);
rank = para(2);
maxIte = para(3);
gamma = para(4);
lambda = para(5);
ite = 0;

U0 = rand(m, rank);
V0 = rand(n, rank);

Lu_plus = (abs(Lu) + Lu) / 2;
Lu_minus = (abs(Lu) - Lu) / 2;

Lv_plus = (abs(Lv) + Lv) / 2;
Lv_minus = (abs(Lv) - Lv) / 2;

while ite <maxIte 
    %[RMSE, MAE] = get_diff(test, U0', V0');
    %fprintf('Ite = %d, RMSE = %0.4f, MAE = %0.4f, time = %0.4f\n', ite, RMSE, MAE, etime(clock, t0));
    %[MAP, MPR, HLU, AUC, avgF, avgP, avgR] = get_diff(test, U0, V0, para);
    %fprintf('Ite = %d, MAP = %0.4f, MPR = %0.4f, HLU = %0.4f, AUC = %0.4f, avgF  = %0.4f, avgP = %0.4f, avgR = %0.4f, time = %0.4f\n', ite, MAP, MPR, HLU, AUC, avgF, avgP, avgR, etime(clock, t0));
    
    U0 = updateU(R, W, P, Lu_plus, Lu_minus, U0, V0, alpha, gamma);
    V0 = updateU(R',W', P', Lv_plus, Lv_minus, V0, U0, alpha, lambda);
    
    ite = ite + 1;
end

U = U0;
V = V0;

end

function [U1] = updateU(R, W, P, Lu_plus, Lu_minus, U0, V, lambda, gamma)
%this updateU function does not use fast equivalence introduced in
%wiZAN_dual paper. This is to control (W)eight and im(P)utation values by
%chem_chem_zinc and protein_protein_zinc_blast indexes.
U1 = U0 .* sqrt( ((W.*(R+P))*V + gamma .* Lu_minus * U0) ./ ( (W.*(U0*V'))*V + gamma.* Lu_plus * U0 + lambda * U0) );
U1(isnan(U1)) = 0;

end
