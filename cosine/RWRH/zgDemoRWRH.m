function zgDemoRWRH(trainCsv, testCsv, dataset)
addpath('/scratch/hansaim.lim/REMAP/benchmark/NTMaxTC/');
addpath('/scratch/hansaim.lim/REMAP/benchmark/NTNL/');
outpath='/scratch/hansaim.lim/wiZAN/cosine/RWRH/aucpr/';
outfile=[outpath dataset '_AUCPR.txt'];

tic;
load('/scratch/hansaim.lim/REMAP/benchmark/sim/chem_chem_zinc.mat');
load('/scratch/hansaim.lim/REMAP/benchmark/sim/prot_prot_zinc.mat');
numChem=size(chem_chem_zinc,1);
numProt=size(prot_prot_zinc,1);
tline=csvread(trainCsv);
Train=sparse(tline(:,1),tline(:,2),1,numChem,numProt);
tsline=csvread(testCsv);
Test=sparse(tsline(:,1),tsline(:,2),1,numChem,numProt);

para = [0.7 0.5 1e-9 100]; %[restart_prob, weight, tolerance, maxIter]
result = zeros(numChem,numProt);

%build graph from similarity matrices and training data
Graph = zeros(numChem+numProt,numChem+numProt);
Graph(1:numChem,1:numChem) = chem_chem_zinc;
Graph(numChem+1:numChem+numProt,numChem+1:numChem+numProt) = prot_prot_zinc;
Graph(1:numChem,numChem+1:numChem+numProt)=Train;
Graph(numChem+1:numChem+numProt,1:numChem)=Train';
%graph built

[nR,nC]=size(Graph); %numRows, numCols
if nR ~= nC
    msg='Graph must be square matrix';
    error(msg);
end

MZ = Graph(1:numChem,1:numChem);
MG = Graph(numChem+1:numChem+numProt,numChem+1:numChem+numProt);
MZG= Graph(1:numChem,numChem+1:numChem+numProt);
MGZ= Graph(numChem+1:numChem+numProt,1:numChem);
Bw1= sum(MZG,2); %sum by Row
Bx = sum(MZ,2); %sum by Row
M  = bsxfun(@rdivide,MZ,Bx); %normalize each row

Bw1(Bw1>0) = 1;
Bw1 = Bw1*(1-para(2));
MZ = bsxfun(@times,M,Bw1);
By = sum(MG,2); %sum by Row
MG = bsxfun(@rdivide,MG,By);
Bw0= sum(MGZ,2);
Bw0(Bw0>0) = 1;
Bw0 = Bw0*(1-para(2));
MG = MG*Bw0;
Bw = sum(MZG,2);
Bz = sum(MGZ,1);
MZG= bsxfun(@rdivide,(para(2)*MZG),Bw);
MGZ= bsxfun(@rdivide,(para(2)*MGZ),Bz);
MZ(1:numChem,1:numChem) = M(1:numChem,1:numChem);
Graph(1:numChem,1:numChem) = MZ;
Graph(numChem+1:numChem+numProt,numChem+1:numChem+numProt) = repmat(MG,1,numProt);
Graph(1:numChem,numChem+1:numChem+numProt) = MZG;
Graph(numChem+1:numChem+numProt,1:numChem) = MGZ;
Graph(isnan(Graph)) = 0;

for i=1:numChem
   restart_vector = zeros(nR,1);
   restart_vector(i,1) = 1;
   r = IterRWR(Graph, para(1), restart_vector, para(4), para(3));
   
   gs = r(numChem+1:numChem+numProt);
   result(i,:) = gs;
end
test_result=TPRbyRowRank(FindTrues(result, Test), 350);
[x,y,tp,auc]=perfcurve(Test(:),result(:),1);
[rec,pre,TPr,AUCpr]=perfcurve(Test(:), result(:), 1, 'xCrit', 'reca', 'yCrit', 'prec');
runTime=toc;
text=['Dataset: ' dataset 'AUC=' num2str(auc) ' AUPR=' num2str(AUCpr) ' Time:' num2str(runTime) '\n'];
fileid=fopen(outfile, 'a+');
fprintf(fileid, text, 'char');
fprintf(fileid,'%7s\t%7s\n','Rank','TPR');
dlmwrite(outfile, test_result, '-append', 'delimiter', '\t', 'precision', 7);
fclose(fileid);
end

function r=IterRWR(Graph, reProb, prefer_vec, maxIter, tolerance)
[nR,nC]=size(Graph);
if nR ~= nC
    msg='Graph is not a square matrix. Check your matrices';
    error(msg);
end
r=prefer_vec;
realIter=maxIter;
ite=0;
while ite <maxIter
    old_r = r;
    r = (1-reProb)*Graph*r + reProb*prefer_vec;
    diff = sum(abs(old_r-r));
    if diff < tolerance
        realIter=ite;
        break
    end
    ite = ite+1;
end
end

function [RCRS] = FindTrues(A, TM) %A: predicted score matrix, TM: Test Matrix (either true positive or true negative)
%returns RCRS=[RowIndex, ColumnIndex, RowRank, PredScore]
[i j k]=find(TM);
RCRS = zeros(length(i),4);
for n=1:length(i)
  row=i(n);
  col=j(n);
  predscore = A(row, col);
  drank = sum(A(row,:) >= A(row, col));	%dense rank on the row
  RCRS(n,:) = [row, col, drank, predscore];
end 
end

function [TPRbyRow] = TPRbyRowRank(RCRS, maxRank)
%to get TPR by cutoff Row rank
%output contains vector of dimension (maxRank, 2)
%Each row of RCRS contains one true positive association: [RowIndex,ColumnIndex,RowRank,PredScore]
%RCRS input must contain row rank on the 3rd column
%RCRS input must contain ONLY test pairs, since the number of lines will be used as condition positive
[cp, col] = size(RCRS); %cp: condition positive, col: number of columns
TPRbyRow = zeros(maxRank, 2);
row_rank = RCRS(:, 3);    %Rank by Rows on 3rd column
rank = 1;       %cutoff rank
while rank <= maxRank
        tpr = sum(row_rank<=rank)/cp;
        TPRbyRow(rank, :) = [rank, tpr];        %update result row
        rank = rank + 1;
end
end
