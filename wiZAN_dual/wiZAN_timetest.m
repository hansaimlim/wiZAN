function wiZAN_timetest(iter, rank, rowsize, colsize)
%iter: number of tests
%rank: factorization low-rank (lower than both row and colsize
%rowsize: number of chemicals
%colsize: number of proteins
tic;
%fixed parameter as of 5/27/2015
para = [0.1, 0.1, 0.01, rank, 400, 0.75, 0.1]; % para: lambda, squared global weight, r, rank, maxIte, gamma, lambda

%chems and prots matrices from chem-chem and prot-prot files
load /scratch/hansaim.lim/wiZAN/ZINC_ChEMBL_DrugBank/chem_chem/chem_chem_ZCD;
load /scratch/hansaim.lim/wiZAN/ZINC_ChEMBL_DrugBank/prot_prot/prot_prot_ZCD;
load /scratch/hansaim.lim/wiZAN/ZINC_ChEMBL_DrugBank/prot_prot/chem_prot_ZCD;
%get number of chemical and protein
m=size(chem_chem_ZCD, 1);	%number of unique chemicals
n=size(prot_prot_ZCD, 1);	%number of unique proteins
%convert csv to matrix
%prots = ceil(prots);
chem_chem_ZCD = chem_chem_ZCD + chem_chem_ZCD';
%prots = prots + prots';
chem_rows=datasample(chem_chem_ZCD,rowsize,1);
chems=datasample(chem_rows,rowsize,2);	%chemical-chemical matrix is square
prot_cols=datasample(prot_prot_ZCD,colsize,1);
prots=datasample(prot_cols,colsize,2);	%protein-protein matrix is square

chemprot_row=datasample(chem_prot_ZCD,rowsize,1);	%contains rowsize-chemicals
chemprots=datasample(chem_prot_row,colsize,2);	%contains colsize-proteins

summ = sum(chems,2); %sum by rows
Dm = spdiags(summ,0,m,m);
Lu = Dm - chems;

sumn = sum(prots,2); %sum by rows
Dn = spdiags(sumn,0,n,n);
Lv = Dn - prots;

[U, V] = updateUV(chemprots, Lu, Lv, para);
P=U*V';	%P is the prediction score matrix
toc
clear
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
