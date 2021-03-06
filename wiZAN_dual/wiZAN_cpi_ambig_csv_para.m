function wiZAN_cpi_ambig_csv_para(true_positive_csv, true_negative_csv, ambiguous_csv, outfile_prefix, para)
%parameter input needed
%para = [0.1, 0.1, 0.01, 300, 100, 0.75, 0.1]; % para: lambda, squared global weight, r, rank, maxIte, gamma, lambda

%chem_chem_zinc and protein_protein_zinc_blast matrices from chem-chem and prot-prot files
load /scratch/hansaim.lim/wiZAN/ZINC_ChEMBL/chem_chem/chem_chem_zinc_chembl;	%loads chem_chem_zc
load /scratch/hansaim.lim/wiZAN/ZINC_ChEMBL/prot_prot/prot_prot_zinc_chembl;	%loads prot_prot_zc
%get number of chemical and protein
temp_c=size(chem_chem_zc);
temp_p=size(prot_prot_zc);
m = temp_c(1);
n = temp_p(1);
%convert csv to matrix
tp_line = csvread(true_positive_csv);
TP = sparse(tp_line(:,1), tp_line(:,2), 1, m, n);      %True Positive activities (i.e. IC50 <= 10uM); ambiguous are removed
tn_line = csvread(true_negative_csv);
TN = sparse(tn_line(:,1), tn_line(:,2), 1, m, n);      %True Negative activities (i.e. IC50 > 10uM); ambiguous are removed
am_line = csvread(ambiguous_csv);
AM = sparse(am_line(:,1), am_line(:,2), 1, m, n);      %AMbiguous activities (multiple tests in both TP and TN range)
clear tp_line;
clear tn_line;
clear am_line;
%protein_protein_zinc_blast = ceil(protein_protein_zinc_blast);
%chem_chem_zc = chem_chem_zc + chem_chem_zc';
%protein_protein_zinc_blast = protein_protein_zinc_blast + protein_protein_zinc_blast';

summ = sum(chem_chem_zc,2); %sum by rows
Dm = spdiags(summ,0,m,m);
Lu = Dm - chem_chem_zc;

sumn = sum(prot_prot_zc,2); %sum by rows
Dn = spdiags(sumn,0,n,n);
Lv = Dn - prot_prot_zc;

[U, V] = updateUV(TP, Lu, Lv, para);
Pred = U*V';	%Predicted score matrix
%known = TP+TN+AM;	%known interactions (including ambiguous ones)
%Pred_unknown = Pred.*(~known);	%Predicted scores for unknown pairs only
%Pred_ambig = Pred.*AM;	%Predicted scores for ambiguous pairs only
%rstat = getRowStat(Pred);
%gstat = getGlobStat(Pred);	%memory heavy
%gstat_un = getGlobStat(Pred_unknown);	%global stat for unknown pairs only
%gstat_am = getGlobStat(Pred_ambig);
%top35_unknown = getTopNbyRow(Pred_unknown, 35);
%rcrsTP = FindTrues(Pred, TP);	%where are True Positives ranked
%rcrsTN = FindTrues(Pred, TN);	%where are True Negatives ranked
%rcrsAM = FindTrues(Pred, AM);
%set output filenames
lrank_u = [outfile_prefix '_lowrank_U.csv'];
lrank_v = [outfile_prefix '_lowrank_V.csv'];
%rstat_filename = [outfile_prefix '_rowstat.csv'];
gstat_filename = [outfile_prefix '_globstat_all.csv'];
%gstat_un_filename = [outfile_prefix '_globstat_unknown.csv'];
%gstat_am_filename = [outfile_prefix '_globstat_ambig.csv'];

%top35_unknown_filename = [outfile_prefix '_top35_unknown.csv'];
%rcrsTP_filename = [outfile_prefix '_TruePositive_rank_score.csv'];
%rcrsTN_filename = [outfile_prefix '_TrueNegative_rank_score.csv'];
%rcrsAM_filename = [outfile_prefix '_AMbiguous_rank_score.csv'];
%set output filenames
csvwrite(lrank_u, U);
csvwrite(lrank_v, V);	%NOT transposed
%csvwrite(rstat_filename, rstat);
csvwrite(gstat_filename, gstat);
%csvwrite(gstat_un_filename, gstat_un);
%csvwrite(gstat_am_filename, gstat_am);
%csvwrite(top35_unknown_filename, top35_unknown);
%csvwrite(rcrTP_filename, rcrTP);
%csvwrite(rcrTN_filename, rcrTN);
%csvwrite(rcrAM_filename, rcrAM);

clear;
end

%function [row_stat] = getRowStat(A) %A is the predicted score matrix; A = U*V'
%si = size(A);
%row_stat = zeros(si(1), 5);
%row_stat = [(1:si(1))', max(A')', min(A')', mean(A')', std(A')', median(A')'];
%end

function [global_stat] = getGlobStat(A) %A is the predicted score matrix
%stats for non-zero values only
[row, col, val] = find(A);
[h, p, adstat, cv] = adtest(val);
global_stat = [max(val), min(val), mean(val), std(val), median(val), h]; %non-normal distribution if h=1
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

function [topN_glob] = getTopNGlobal(A, N)	%A: input matrix, N: how many?
si = size(A);
topN_glob = zeros(N,3);	%[row, col, val]
count = 1;
[sortval, sorti] = sort(A(:), 'descend');
for idx = [sorti(1:N)]
 [row, col] = ind2sub(si, idx);
 val = A(row,col);
 topN_glob(count,:) = [row, col, val];
 count = count + 1;
end
end

function [row_col_rank_score] = FindTrues(A, TM) %A: predicted score matrix, TM: True Matrix (either true positive or true negative)
si = (size(A));
row_col_rank_score = zeros(sum(TM(:)),4);
rcrs_count = 1;
for row = 1:si(1)
 true_index = find(TM(row,:));
 for ind = 1:length(true_index)
  col = true_index(ind);
  drank = sum(A(row,:) >= A(row, col));
  predscore = A(row, col);
  row_col_rank_score(rcr_count,:) = [row, col, drank, predscore];
  rcrs_count = rcrs_count + 1;
 end
end

end

function [U, V] = updateUV(R, Lu, Lv, para)
% para: lambda, r, T, rank, maxIte, ite_of_bisection method, topN
[m, n] = size(R);
alpha = para(1);
w = para(2);
r = para(3);
rank = para(4);
maxIte = para(5);
gamma = para(6);
lambda = para(7);
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
    U0 = updateU(R, UVT, w, r, Lu_plus, Lu_minus, U0, V0, alpha, gamma);
    V0 = updateU(R', UVT',w, r, Lv_plus, Lv_minus, V0, U0, alpha, lambda);
    
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

function [U1] = updateU(R, UVT, w, p, Lu_plus, Lu_minus, U0, V, lambda, gamma)

[m,n]=size(R);
U1 = U0 .* sqrt( ((1-w*p)*R*V + ones(m,1)*p*((w*ones(1,n))*V) + gamma .* Lu_minus * U0) ./ ((1-w)*UVT*V + w * (U0*(V'*V))  + gamma.* Lu_plus * U0 + lambda * U0) );

U1(isnan(U1)) = 0;

end
