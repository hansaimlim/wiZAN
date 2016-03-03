function wiZAN_cpi_controlWP
true_positive_csv='/scratch/hansaim.lim/wiZAN/ZINC_ChEMBL/chem_prot/zc_active.csv';
true_negative_csv='/scratch/hansaim.lim/wiZAN/ZINC_ChEMBL/chem_prot/zc_inactive.csv';
%ambiguous_csv='/scratch/hansaim.lim/wiZAN/ZINC_ChEMBL/chem_prot/zc_ambiguous.csv';
outfile_dir='/scratch/hansaim.lim/wiZAN_on_zinc_chembl/WPcontrol/';
outfile_prefix='zc_controlWP';

%fixed parameter as of 5/27/2015
%parameter(5) updated to 400 (iterations)
%W and P matrices are controlled by true negative pairs
%W=1, P=0, R=0 if the pair is inactive (true negative)
%otherwise, w=0.1, p=0.01
para = [0.1, 300, 400, 0.75, 0.1]; % para: lambda, rank, maxIte, gamma, lambda

%chem_chem_zinc and protein_protein_zinc_blast matrices from chem-chem and prot-prot files
load /scratch/hansaim.lim/wiZAN/ZINC_ChEMBL/chem_chem/chem_chem_zinc_chembl;	%loads chem_chem_zc
load /scratch/hansaim.lim/wiZAN/ZINC_ChEMBL/prot_prot/prot_prot_zinc_chembl;	%loads prot_prot_zc
%get number of chemical and protein
temp_c=size(chem_chem_zc);
temp_p=size(prot_prot_zc);
m = temp_c(1);
n = temp_p(1);
%convert csv to matrix
%protein_protein_zinc_blast = ceil(protein_protein_zinc_blast);
%chem_chem_zc = chem_chem_zc + chem_chem_zc';
%protein_protein_zinc_blast = protein_protein_zinc_blast + protein_protein_zinc_blast';

summ = sum(chem_chem_zc,2); %sum by rows
Dm = spdiags(summ,0,m,m);
Lu = Dm - chem_chem_zc;

sumn = sum(prot_prot_zc,2); %sum by rows
Dn = spdiags(sumn,0,n,n);
Lv = Dn - prot_prot_zc;

tpline=csvread(true_positive_csv);
tnline=csvread(true_negative_csv);
TP=sparse(tpline(:,1), tpline(:,2), 1, m, n);
TN=sparse(tnline(:,1), tnline(:,2), 1, m, n);
W=(ones(m,n).*(~(TP+TN))).*0.1 + TN;	%0.1 if unobserved, 1 if unassociated, 0 otherwise
P=(ones(m,n).*(~(TP+TN))).*0.01;	%0.01 if unobserved, 0 otherwise
 tic;
 [U, V] = updateUV(TP, Lu, Lv, para, W, P);
 toc
% Pred = U*V';	%Predicted score matrix
% rcrs = FindTrues(Pred, TS);	%get the ranks for Test pairs

lrankU=[ outfile_dir outfile_prefix '_lowrank_U.csv'];
lrankV=[ outfile_dir outfile_prefix '_lowrank_V.csv'];
csvwrite(lrankU, U);
csvwrite(lrankV, V);
clear;
end

function [row_stat] = getRowStat(A) %A is the predicted score matrix; A = U*V'
si = size(A);
row_stat = zeros(si(1), 5);
row_stat = [(1:si(1))', max(A')', min(A')', mean(A')', std(A')'];
end

function [global_stat] = getGlobStat(A) %A is the predicted score matrix
global_stat = [max(A(:)), min(A(:)), mean(A(:)), std(A(:))];
end

function [topN_byRow] = getTopNbyRow(A, N) %A: predicted score matrix, N: how many you want to see from the top
%topN_byRow contains column indices from top rank(by score) to Nth rank
%row indices are equal to the line number
si = (size(A));
topN_byRow = zeros(si(1), N);
for K = 1:si(1)
 [sortval, sorti] = sort(A(K,:), 'descend');
 topN_byRow(K,:) = sorti(1:N);
end
end

function [U, V] = updateUV(R, Lu, Lv, para, W, P)
% para: lambda, rank, maxIte, ite_of_bisection method, importance on similarity matrices
[m, n] = size(R);
alpha = para(1);
rank = para(2);
maxIte = para(3);
gamma = para(4);
lambda = para(5);
ite = 0;

%IU = ones(m,n) - R;
%W = IU * w + R;
%IU = sparse(IU) * r;

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
    

    UVT = get_UVT(R, U0, V0);
    U0 = updateU(R, W, P, Lu_plus, Lu_minus, U0, V0, alpha, gamma);
    V0 = updateU(R', W', P', Lv_plus, Lv_minus, V0, U0, alpha, lambda);
    
    ite = ite + 1;
end

U = U0;
V = V0;

end

function [UVT] = get_UVT(R, U, V)

[m, n] = size(R);

[I, J] = find(R);
iSize = size(I, 1);
K = ones(iSize,1);
count = 0;
for j = 1:n
    id = find(R(:,j));
    len = length(id);
    K(count+1:count+len) = U(id, :) * V(j,:)';
    count = count + len;
end

UVT = sparse(I, J, K, m, n);

end

function [U1] = updateU(R, W, P, Lu_plus, Lu_minus, U0, V, lambda, gamma)

U1 = U0 .* sqrt( ((W.*(R+P))*V + gamma .* Lu_minus * U0) ./ ( (W.*(U0*V'))*V + gamma.* Lu_plus * U0 + lambda * U0) );

U1(isnan(U1)) = 0;

end

function [U1] = updateU_lowrank(R, UVT, w, p, Lu_plus, Lu_minus, U0, V, lambda, gamma)

[m,n]=size(R);
U1 = U0 .* sqrt( ((1-w*p)*R*V + ones(m,1)*p*((w*ones(1,n))*V) + gamma .* Lu_minus * U0) ./ ((1-w)*UVT*V + w * (U0*(V'*V))  + gamma.* Lu_plus * U0 + lambda * U0) );

U1(isnan(U1)) = 0;

end
