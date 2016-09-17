function REMAP_chem_dis_lin()
%[1] Hansaim Lim, Aleksandar Poleksic, Hanghang Tong, Yuan Yao, Di He, Luke Zhuang, Patrick Meng, and Lei Xie, "Large-scale Off-target Identification Using Fast and Accurate Dual Regularized One-Class Collaborative Filtering and Its Application to Drug Repurposing" , under review

%10fold cross validation for ZINC dataset
%input_dir='../benchmark/cv10/';
%chem_chem_zinc and protein_protein_zinc_blast matrices from chem-chem and prot-prot files
%load ../benchmark/sim/chem_chem_zinc;
%load ../benchmark/sim/protein_protein_zinc_blast;
%get number of chemical and protein
load('C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\data\similarity\chem_chem_sim.mat');
load('C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\data\similarity\prot_prot_sim.mat');
load('C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\data\similarity\dis_dis_sim_lin_norm.mat');
load('C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\data\similarity\dis_dis_sim_path_norm.mat');
load('C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\data\similarity\dis_dis_sim_resnik_norm.mat');
load('C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\data\similarity\dis_dis_sim_vector_norm.mat');
%load /scratch/hansaim.lim/wiZAN/fascinate/data/chem_gene/chem_gene
%load /scratch/hansaim.lim/wiZAN/fascinate/data/chem_dis/chem_dis
%load /scratch/hansaim.lim/wiZAN/fascinate/data/gene_dis/gene_dis
outfile='C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\REMAP\dis_sim\TPR20_lin.txt';
input_dir='C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\data\chem_dis\10fold\';
chem_chem_mat=chem_chem_sim;
prot_prot_mat=dis_dis_sim_lin_norm;
m=size(chem_chem_mat, 1);
n=size(prot_prot_mat, 1);

summ = sum(chem_chem_mat,2); %sum by rows
Dm = spdiags(summ,0,m,m);
Lu = Dm - chem_chem_mat;

sumn = sum(prot_prot_mat,2); %sum by rows
Dn = spdiags(sumn,0,n,n);
Lv = Dn - prot_prot_mat;

disp(['Quick REMAP test on chem-dis dataset'])
rank=75;
iter=300;
para = [0.1, 0.1, 0.01, rank, iter, 0.75, 0.1]; % para: lambda, squared global weight, r, rank, maxIte, gamma, lambda

TPR20=zeros(1,10);
for k=1:10
 trainfile=[input_dir 'train' num2str(k) '.csv'];
 testfile =[input_dir 'test' num2str(k) '.csv'];
 trline=csvread(trainfile);
 tsline=csvread(testfile);
 TR=sparse(trline(:,1), trline(:,2), 1, m, n);
 TS=sparse(tsline(:,1), tsline(:,2), 1, m, n);
 [U, V] = updateUV(TR, Lu, Lv, para);
 test_result = TPRbyRowRank(FindTrues(U*V', TS), 20);   %max cutoff rank 20
 TPR20(1,k)=test_result(20,2);
 clear U V TR TS trline tsline test_result;
end
fid=fopen(outfile,'a+');
fprintf(fid,['rank=' num2str(rank) ' iter=' num2str(iter) ' dis_sim=lin_normalized ' 'TPR20\n']);
dlmwrite(outfile,TPR20,'-append');
TPR20
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
