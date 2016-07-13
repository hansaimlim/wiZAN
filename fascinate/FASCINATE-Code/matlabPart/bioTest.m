% the experiment is based on BIO data set
clear
%load construct chemical layer
load chemSim
maxi = max(I);
maxj = max(J);
nChem = max(maxi,maxj);
Achem = sparse(I,J,weight,nChem,nChem);
%construct disease layer
load diseaseSim
maxi = max(I);
maxj = max(J);
nDis = max(maxi,maxj);
Adis = sparse(I,J,weight,nDis,nDis);
%construct gene layer
load geneSim
maxi = max(I);
maxj = max(J);
nGene = max(maxi,maxj);
Agene = sparse(I,J,ones(length(I),1),nGene,nGene);
%load and construct cross-layer dependencies
load chem_disDep
chem_dis = sparse(chem,disease,ones(length(chem),1),nChem,nDis);
load chem_geneDep
chem_gene = sparse(chem, gene,ones(length(chem),1),nChem,nGene);
load gene_disDep
gene_dis = sparse(gene,disease,ones(length(gene),1),nGene,nDis);

% make similarity matrix a symmetric matrix
Achem = Achem + Achem';
Adis = Adis + Adis';
Agene =  Agene + Agene';

%make each layer and cross-layer dependency matrices to be unweighted graphs(optional)
Agene(Agene>1) = 1;
chem_dis(chem_dis>1) = 1;
chem_gene(chem_gene>1) = 1;
gene_dis(gene_dis>1) = 1;



n = 3; %number of graphs
%initialize all matrices
G{1}.A = Achem;
G{2}.A = Agene;
G{3}.A = Adis;
%define the graph-level dependency relationship
D = zeros(n,n);
D(1,2) = 1;
D(1,3) = 1;
D(2,3) = 1;
%make  hashmap D_new and DU
%construct 'hashmap' with (row,col) and number
[row,col] = find(D);
links = length(row);
D_new = sparse(row,col,1:links,size(D,1),size(D,2));
D_new = D_new+D_new';

%define node-level dependency relationship
DU{1}.D = chem_gene;
DU{2}.D = chem_dis;
DU{3}.D = gene_dis;



%%%%%%%%%%%%%%%%
% start inference
%%%%%%%%%%%%%%%%%%%
fprintf('start inference');
[ F ] = mulanImputeF( G,D_new,DU,0.1,0.1);

%construct the inferred results
Infer{1}.D = F{1}.F*F{2}.F';
Infer{2}.D = F{1}.F*F{3}.F';
Infer{3}.D = F{2}.F*F{3}.F';

%filter to get the scores of the missing dependencies only
for i = 1:length(Infer)
    [I,J] = find(DU{i}.D);
    mask = sparse(I,J,ones(length(I),1),size(DU{i}.D,1),size(DU{i}.D,2));
    Score{i}.D = Infer{i}.D - mask.*Infer{i}.D;
end
