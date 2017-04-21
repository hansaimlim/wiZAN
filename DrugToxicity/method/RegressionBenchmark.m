function RegressionBenchmark()
%Benchmark Regression code on Affymetrix Liver gene expression dataset
%use 10-fold Cross Validation
%varying number of splits and learning rates
%Generates plots of MSE vs. NumTrees
addpath('D:/OneDrive - Cuny GradCenter/DrugMatrix Manual/AffyLiver/');
load('AffyLiver.mat');
load('AffyLiverMeanResponse.mat');
InputMat=AffyLiver;
Observations=AffyLiverMeanResponse(:,1);
k=10; %kFold Cross validation
% InputMat=AffyLiver(1:20,:);
% Observations=AffyLiverMeanResponse(1:20,1);

maxNumCompThreads(3); %Max. number of cores to use

n = size(AffyLiver,1);
m = floor(log2(n - 1));
learnRate = [0.1 0.25 0.5 0.75 1];
numLR = numel(learnRate);
maxNumSplits = 2.^(0:m/2);
numMNS = numel(maxNumSplits);
numTrees = 150;
% Mdl = cell(numMNS,numLR);
c=cvpartition(size(InputMat,1),'kFold',k);
predDeep=zeros(size(InputMat,1),1);
predStump=zeros(size(InputMat,1),1);
for kF=1:k
    idxTrain=training(c,kF);
    idxTest=test(c,kF);
    TrainInput=InputMat(idxTrain,:);
    TestInput=InputMat(idxTest,:);
    TrainObs=Observations(idxTrain,1);
    TestObs=Observations(idxTest,1);

    %Baseline Regression
    MdlDeep = fitrtree(TrainInput,TrainObs,'MergeLeaves','off',...
    'MinParentSize',1,'Surrogate','on');
    MdlStump = fitrtree(TrainInput,TrainObs,'MaxNumSplits',1,...
    'Surrogate','on');
    predDeep(idxTest)=predict(MdlDeep,TestInput);
    predStump(idxTest)=predict(MdlStump,TestInput);
end
[CorDeep,PCorDeep]=Correlation(predDeep,Observations)
rsqDeep=Rsquared(predDeep,Observations)
[CorStump,PCorStump]=Correlation(predStump,Observations)
rsqStump=Rsquared(predStump,Observations)


for i = 1:numLR %LearnRate
    LR=learnRate(i);
    CORR=zeros(numMNS*numTrees,3);
    PCOR=zeros(numMNS*numTrees,3);
    RSQ=zeros(numMNS*numTrees,3);
    count=1;
    for j = 1:numMNS %MaxNumSplits
        mns=maxNumSplits(j);
        for nt=1:numTrees
            pred=zeros(size(InputMat,1),1);
            for kF=1:k
                idxTrain=training(c,kF);
                idxTest=test(c,kF);
                TrainInput=InputMat(idxTrain,:);
                TestInput=InputMat(idxTest,:);
                TrainObs=Observations(idxTrain,:);
                TestObs=Observations(idxTest,:);
                
                %Ensemble Trees
                t = templateTree('MaxNumSplits',mns,'Surrogate','on');
                Mdl = fitrensemble(TrainInput,TrainObs,'NumLearningCycles',nt,...
                'Learners',t,'LearnRate',LR);
                pred(idxTest)=predict(Mdl,TestInput);
            end
            %correlation coefficient
            [cor,pval_cor]=Correlation(pred,Observations);
            CORR(count,:)=[nt,mns,cor];
            PCOR(count,:)=[nt,mns,pval_cor];

            %R-squared
            Rsq=Rsquared(pred,Observations);
            RSQ(count,:)=[nt,mns,Rsq];
            count=count+1;
        end
    end
    outfile=['./RegressionBenchmark_LearnRate' num2str(LR) '_correlation.txt'];
    fid = fopen(outfile,'wt');
    fprintf(fid, '%s\t%s\t%s\n', 'numTrees','maxNumSplits','Correlation');  % header
    fclose(fid);
    dlmwrite(outfile,CORR,'-append','delimiter','\t','precision',10);
    
    outfile=['./RegressionBenchmark_LearnRate' num2str(LR) '_correlationPvalue.txt'];
    fid = fopen(outfile,'wt');
    fprintf(fid, '%s\t%s\t%s\n', 'numTrees','maxNumSplits','Correlation');  % header
    fclose(fid);
    dlmwrite(outfile,PCOR,'-append','delimiter','\t','precision',10); 
    
    outfile=['./RegressionBenchmark_LearnRate' num2str(LR) '_Rsquared.txt'];
    fid = fopen(outfile,'wt');
    fprintf(fid, '%s\t%s\t%s\n', 'numTrees','maxNumSplits','Correlation');  % header
    fclose(fid);
    dlmwrite(outfile,RSQ,'-append','delimiter','\t','precision',10);
    
end

end