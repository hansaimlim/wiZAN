function mulanImputeF_10cv(input_dir,option,alpha,beta)
%MULANIMPUTE Summary of this function goes here
%   the function takes a multi-layered network as input and output the low
%   rank approximation for each layer
% \sum||W.(D(i,j)-Fi'Fj))||^2 +\sum \beta tr(Fi'(Ti-Ai)Fi)+ \sum beta||Fi||^2
% INPUT: G: adjacency matrix for each layers; D_map: dependency map matrix;
% DM: dependency matrice; beta:coefficient of ||Fi||; alpha:coefficients
% of user homopily;
% OUTPUT: set of low-rank approximations

% input_dir is the path to directory where 10-folded .csv files are stored
% input file names must be trainX.csv and testX.csv where X is the integer number from 1 to 10
% option is an integer tag (1,2,3)
% option=1 : cross-validate chemical-gene
% option=2 : cross-validate chemical-disease
% option=3 : cross-validate gene-disease

load /scratch/hansaim.lim/wizan/wiZAN/fascinate/data/similarity/chem_chem_sim
load /scratch/hansaim.lim/wizan/wiZAN/fascinate/data/similarity/prot_prot_sim
load /scratch/hansaim.lim/wizan/wiZAN/fascinate/data/similarity/dis_dis_sim_blank
load /scratch/hansaim.lim/wizan/wiZAN/fascinate/data/chem_gene/chem_gene
load /scratch/hansaim.lim/wizan/wiZAN/fascinate/data/chem_dis/chem_dis
load /scratch/hansaim.lim/wizan/wiZAN/fascinate/data/gene_dis/gene_dis
G{1}.A=chem_chem_sim;
G{2}.A=prot_prot_sim;
G{3}.A=dis_dis_sim;
D_new=sparse([1,1,2],[2,3,3],[1,2,3],3,3);
D_new=D_new+D_new';
topnum=0;%the number of top1% prediction varies according to the option
testcase='';
row=0;
col=0; %row and column sizes for train and test file. updated according to the option
if option==1
    %option for chem-gene
    topnum=floor((min(size(chem_chem_sim,1),size(prot_prot_sim,1)))/100);
    testcase='chem-gene';
    row=size(chem_chem_sim,1);
    col=size(prot_prot_sim,1);
elseif option==2
    %option for chem-dis
    topnum=floor((min(size(chem_chem_sim,1),size(dis_dis_sim,1)))/100);
    testcase='chem-dis';
    row=size(chem_chem_sim,1);
    col=size(dis_dis_sim,1);
elseif option==3
    %option for gene-dis
    topnum=floor((min(size(prot_prot_sim,1),size(dis_dis_sim,1)))/100);
    testcase='gene-dis';
    row=size(prot_prot_sim,1);
    col=size(dis_dis_sim,1);
else
    %option not properly chosen
    %print error message and terminate
    msg='Please set an option. 1 for chem-gene, 2 for chem-dis, or 3 for gene-dis validation.\n';
    error(msg);
end
if topnum<20
    %top 1% is too small
    topnum=20;
end
%define dependencies for all options
DU{1}.D=chem_gene;
DU{2}.D=chem_dis;
DU{3}.D=gene_dis;
%chosen option will be overwritten later (e.g. chem-dis by training set if option==2)

paras=[0.1,100,400]; % weight, rank, maxIte
rank=paras(2);
tic;
%the same initial low-rank matrices for the same rank
for i = 1:length(G)
    m = size(G{i}.A,1);
    F0{i}.F =  (rand(m, rank))/sqrt(rank);
end

TPR_sum=0;
for k=1:10
    trainfile=[input_dir 'train' num2str(k) '.csv'];
    testfile=[input_dir 'test' num2str(k) '.csv'];
    trline=csvread(trainfile);
    tsline=csvread(testfile);
    TR=sparse(trline(:,1),trline(:,2),1,row,col);
    TS=sparse(tsline(:,1),tsline(:,2),1,row,col);

    if option==1
	%option for chem-gene
	%train matrix instead of DU{1}.D
	DU{1}.D=TR;
    elseif option==2
	%option for chem-dis
	%train matrix instead of DU{2}.D
	DU{2}.D=TR;
    elseif option==3
	%option for gene-dis
	%train matrix instead of DU{3}.D
	DU{3}.D=TR;
    else
	%option not properly chosen
	%print error message and terminate
	msg='Please set an option. 1 for chem-gene, 2 for chem-dis, or 3 for gene-dis validation.\n';
	error(msg);
    end
    [F] = updateF(G,D_new,DU,beta,alpha,F0,paras);
    P=zeros(2,2); %a temporary matrix for prediction results. updated according to the given option
    if option==1
	%option for chem-gene
	P=F{1}.F*F{2}.F';
    elseif option==2
	%option for chem-dis
	P=F{1}.F*F{3}.F';
    elseif option==3
	%option for gene-dis
	P=F{2}.F*F{3}.F';
    else
	%option not properly chosen
	%print error message and terminate
	msg='Please set an option. 1 for chem-gene, 2 for chem-dis, or 3 for gene-dis validation.\n';
	error(msg);
    end
    test_result=TPRbyRowRank(FindTrues(P,TS),topnum);%max cutoff rank is topnum (defined according to option)
    TPR_sum=TPR_sum+test_result(topnum,2);
    clear trainfile testfile TR TS trline tsline test_result;
end
TPR_avg=(TPR_sum/10); %Avg. TPR value for 10-fold CV
disp(['10CV TPR' num2str(topnum) ': ' num2str(v) 'at rank=' rank ', iter=' num2str(paras(3)) ' option=' testcase ' alpha=' num2str(alpha) ' beta=' num2str(beta)])
toc

end


function [F] = updateF(G,D_map,DM,beta,alpha,F0,paras)
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
