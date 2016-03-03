function wiZAN_singletest(chem_chem_sim, prot_prot_sim, chem_prot, outfile_dir, outfile_prefix)
tic;	%measuring the running time
%fixed parameter as of 5/27/2015
%parameter(5) updated to 400 (iterations)
para = [0.1, 0.1, 0.01, 300, 400, 0.75, 0.1]; % para: lambda, squared global weight, r, rank, maxIte, gamma, lambda

%chem_chem_zinc and protein_protein_zinc_blast matrices from chem-chem and prot-prot files
%get number of chemical and protein
m=size(chem_chem_sim,1);	%number of chemicals
n=size(prot_prot_sim,1);	%number of proteins
%convert csv to matrix
%protein_protein_zinc_blast = ceil(protein_protein_zinc_blast);
%chem_chem_sim = chem_chem_sim + chem_chem_sim';
%protein_protein_zinc_blast = protein_protein_zinc_blast + protein_protein_zinc_blast';

summ = sum(chem_chem_sim,2); %sum by rows
Dm = spdiags(summ,0,m,m);
Lu = Dm - chem_chem_sim;

sumn = sum(prot_prot_sim,2); %sum by rows
Dn = spdiags(sumn,0,n,n);
Lv = Dn - prot_prot_sim;

[U, V] = updateUV(chem_prot, Lu, Lv, para);
Pred = U*V';	%Predicted score matrix
toc	%First running time-how long the program runs
tic;	%second running time-for writing file output
lrank_u = [outfile_dir outfile_prefix '_lowrank_U.csv'];
lrank_v = [outfile_dir outfile_prefix '_lowrank_V.csv'];
csvwrite(lrank_u, U);
csvwrite(lrank_v, V);
toc	%second running time-for writing file output
clear;
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

