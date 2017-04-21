function OptimizeRandomForest(InputMat,Observations,kFold)
maxNumCompThreads(3); %determine the maximum number of cores to use
k=kFold;
c=cvpartition(size(InputMat,1),'kFold',k);
numTrees=20:20:500;
numPredictors=1:floor((size(InputMat,2)/100)):floor(size(InputMat,2)/10);
% numTrees=1:5; %quick test
% numPredictors=10:10:50; %quick test
%avgMSE=zeros(numel(numTrees),numel(numPredictors));
%stdMSE=zeros(numel(numTrees),numel(numPredictors));
COR=zeros(numel(numTrees),numel(numPredictors));
PCOR=zeros(numel(numTrees),numel(numPredictors));
RSQ=zeros(numel(numTrees),numel(numPredictors));
for nt=1:numel(numTrees)
    numTree=numTrees(nt);
    for np=1:numel(numPredictors)
        numPredictor=numPredictors(np);
        %MSE=zeros(1,k);
        pred=zeros(size(InputMat,1),1);
        for kF=1:k
            idxTrain=training(c,kF);
            idxTest=test(c,kF);
            TrainInput=InputMat(idxTrain,:);
            TestInput=InputMat(idxTest,:);
            TrainObs=Observations(idxTrain,1);
            TestObs=Observations(idxTest,1);

            %Default linear reg.
            Mdl=TreeBagger(numTree,TrainInput,TrainObs,'Method','regression','NumPredictorsToSample',numPredictor);
            pred(idxTest)=predict(Mdl,TestInput);          
            %L=error(Mdl,TestInput,TestObs,'Mode','Ensemble');
            %MSE(1,kF)=L; %MSE from default linear reg.
        end
        %avgMSE(nt,np)=mean(MSE);
        %stdMSE(nt,np)=std(MSE);
        [Corr,PCorr]=Correlation(pred,Observations);
        Rsq=Rsquared(pred,Observations);
        COR(nt,np)=Corr;
        PCOR(nt,np)=PCorr;
        RSQ(nt,np)=Rsq;
    end
end
CORRS=zeros(numel(numTrees)*numel(numPredictors),3);
PCORRS=zeros(numel(numTrees)*numel(numPredictors),3);
RSQS=zeros(numel(numTrees)*numel(numPredictors),3);
count=1;
for nt=1:size(COR,1)
    numT=numTrees(nt);
    for np=1:size(COR,2)
        numP=numPredictors(np);
        cor=COR(nt,np);
        pcor=PCOR(nt,np);
        rsq=RSQ(nt,np);
        CORRS(count,:)=[numT,numP,cor];
        PCORRS(count,:)=[numT,numP,pcor];
        RSQS(count,:)=[numT,numP,rsq];
        count=count+1;
    end
end
clear count nt np numT numP val

fid = fopen('./OptimizeRF_Correlations.txt','wt');
fprintf(fid, '%s\t%s\t%s\n', 'numTrees','numPredictors','Correlation');  % header
fclose(fid);
dlmwrite('./OptimizeRF_Correlations.txt',CORRS,'-append','delimiter','\t','precision',10);

fid = fopen('./OptimizeRF_CorrelationPvalues.txt','wt');
fprintf(fid, '%s\t%s\t%s\n', 'numTrees','numPredictors','Correlation_Pvalue');  % header
fclose(fid);
dlmwrite('./OptimizeRF_CorrelationPvalues.txt',PCORRS,'-append','delimiter','\t','precision',10);

fid = fopen('./OptimizeRF_Rsquared.txt','wt');
fprintf(fid, '%s\t%s\t%s\n', 'numTrees','numPredictors','RSquared');  % header
fclose(fid);
dlmwrite('./OptimizeRF_Rsquared.txt',RSQS,'-append','delimiter','\t','precision',10);


end
