function [CG, CS, GS]=Fascinate_Novartis()
%MULANIMPUTE Summary of this function goes here
%   the function takes a multi-layered network as input and output the low
%   rank approximation for each layer
% \sum||W.(D(i,j)-Fi'Fj))||^2 +\sum \beta tr(Fi'(Ti-Ai)Fi)+ \sum beta||Fi||^2
% INPUT: G: adjacency matrix for each layers; D_map: dependency map matrix;
% DM: dependency matrice; beta:coefficient of ||Fi||; alpha:coefficients
% of user homopily;
% OUTPUT: set of low-rank approximations

% This script is to optimize rank and iteration parameter based on 10-fold cross validation
% optimize performance for chem-gene and chem-dis prediction to find best
% parameters for the two layers
% Using the optimized parameters, get gene-dis prediction

% Do not discard the low rank matrices

%optimal parameters are:
rank=770; iter=110; alpha=0.5; beta=1.0; wt=0.25;

maxNumCompThreads(12); %determine the maximum number of cores to use
kFold=10; %use 10-fold CV to optimize parameters

% Real data
 load('/scratch/hansaim.lim/wiZAN/fascinate/data/ppi_biogrid.mat');
 load('/scratch/hansaim.lim/wiZAN/sider/data/similarity/chem_chem_SIDER.mat');
 load('/scratch/hansaim.lim/wiZAN/sider/data/similarity/se_se_blank.mat');
 load('/scratch/hansaim.lim/wiZAN/sider/data/association/gene_se_Novartis.mat');
 load('/scratch/hansaim.lim/wiZAN/sider/data/association/chem_gene_SIDER.mat');
 load('/scratch/hansaim.lim/wiZAN/sider/data/association/chem_se_SIDER.mat');
 chem_chem=chem_chem_SIDER;
 prot_prot=ppi_biogrid;
 se_se=se_se_blank;
 chem_gene=chem_gene_SIDER;
 chem_se=chem_se_SIDER;
 gene_se=gene_se_Novartis;
% Real data

%dummy data for script testing
%chem_chem=rand(50,50);chem_chem=(chem_chem+chem_chem')./2;
%prot_prot=rand(60,60);prot_prot=(prot_prot+prot_prot')./2;
%se_se=rand(40,40);se_se=(se_se+se_se')./2;
%chem_gene=rand(50,60);chem_gene(chem_gene<0.7)=0;chem_gene(chem_gene>0)=1;
%chem_se=rand(50,40);chem_se(chem_se<0.7)=0;chem_se(chem_se>0)=1;
%gene_se=rand(60,40);gene_se(gene_se<0.7)=0;gene_se(gene_se>0)=1;
%dummy data for script testing

G{1}.A=chem_chem;
G{2}.A=prot_prot;
G{3}.A=se_se;
D_new=sparse([1,1,2],[2,3,3],[1,2,3],3,3);
D_new=D_new+D_new';
%define dependencies for all options
DU{1}.D=chem_gene;
DU{2}.D=chem_se;
DU{3}.D=gene_se;

for ig = 1:length(G)
    m = size(G{ig}.A,1);
    F0{ig}.F =  (rand(m, rank))/sqrt(rank);
end

paras=[wt,rank,iter]; %[weight, rank, iter]
[F] = updateF(G,D_new,DU,alpha,beta,F0,paras);
Chemlow=F{1}.F;
Genelow=F{2}.F;
SElow=F{3}.F;
save('/scratch/hansaim.lim/wiZAN/fascinate/Novartis/Chem_LowRank.mat','Chemlow','-v7.3');
save('/scratch/hansaim.lim/wiZAN/fascinate/Novartis/Gene_LowRank.mat','Genelow','-v7.3');
save('/scratch/hansaim.lim/wiZAN/fascinate/Novartis/SE_LowRank.mat','SElow','-v7.3');
disp('LowRank matrices are saved under /scratch/hansaim.lim/wiZAN/fascinate/Novartis/ \n');
CG=Chemlow*Genelow'; %chem-gene predicted
CS=Chemlow*SElow'; %chem-se predicted
GS=Genelow*SElow'; %gene-se predicted
end


function [F] = updateF(G,D_map,DM,alpha,beta,F0,paras)
[I,J,V] = find(D_map);
weight = paras(1);
maxIter = paras(3);

threshold = 1e-8;

%initialized low rank approximations
for i = 1:length(G)
    F{i}.F = F0{i}.F;
end


maxChange = 1;
iter = 0;

%iteratively update each low-rank matrix
while maxChange > threshold && iter < maxIter
    %update each low-rank approximation matrix
    for i = 1:length(G)     
        i;
        F0 = F{i}.F;
        F{i}.F = updateFi(i,D_map,weight,G,DM,F,beta,alpha);
        changes(i) = max(max(abs(F{i}.F-F0)));
    end
    maxChange  = max(changes);
    iter = iter + 1;
end
iter;
end

function [Fi_new] = updateFi(curLayer,Dmap,weight,G,DM,F,beta,alpha)

Fi = F{curLayer}.F;
[m,n] = size(Fi);
%initialize A2 and B2
A2 = zeros(m,n);
B2 = zeros(m,n);
%find related layers
relLayers = find(Dmap(curLayer,:));
%get summation part of A2 and B2
for i = 1:length(relLayers)
    j = relLayers(i);
    index = Dmap(curLayer,j);
    Dij = DM{index}.D;
    if curLayer > relLayers(i)
        Dij = Dij';
    
    end
    A2 = A2 + Dij*F{j}.F;
    
    Rhat = get_UVT(Dij,Fi,F{j}.F);
    B2 = B2 + (1-weight)*Rhat*F{j}.F+weight*Fi*(F{j}.F'*F{j}.F);
end
%remaining part of A2
A2 = A2 + alpha*G{curLayer}.A*Fi;
%remaining parts of B2
sumA = sum(G{curLayer}.A,2);
DA = spdiags(sumA,0,m,m);
B2 = B2 + beta*Fi + alpha*DA*Fi + 1e-8;

Fi_new = Fi.*sqrt(A2./B2);

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
