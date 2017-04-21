function OptimizeGLMNet(InputMat,Observations,family,kFold)
maxNumCompThreads(3); %determine the maximum number of cores to use
%Usage:
%OptimizaGLMNet(InputMat,Observations,distribution_family,kFold,output_file);
%For example:
%OptimizeGLMNet(AffyLiver,AffyLiverMeanResponse(:,1),'gaussian',10,'./GLMGaussian_AffyMeanResponse.tsv');
%OptimizeGLMNet(AffyLiver,AffyLiverMeanResponse(:,1),'poisson',10,'./GLMPoisson_AffyMeanResponse.tsv');

addpath('C:\Users\hansaimlim\OneDrive - Cuny GradCenter\DrugMatrix Manual\algo\LASSO\glmnet_matlab');
% load('C:\Users\hansaimlim\OneDrive - Cuny GradCenter\DrugMatrix Manual\AffyLiver\AffyLiver.mat');
% load('C:\Users\hansaimlim\OneDrive - Cuny GradCenter\DrugMatrix Manual\AffyLiver\AffyLiverMeanResponse.mat');
% InputMat=AffyLiver;
% Observations=AffyLiverMeanResponse(:,1);

k=kFold;
if isempty(family)
    family = 'gaussian';
end
c=cvpartition(size(InputMat,1),'kFold',k);

alpha=0:0.1:1;
nlambda=1:floor((size(InputMat,2)/100)):floor(size(InputMat,2)/10);
% alpha=0:0.5:1;
% nlambda=50:50:150;
% avgMSE=zeros(numel(alpha),numel(nlambda));
% stdMSE=zeros(numel(alpha),numel(nlambda));
COR=zeros(numel(alpha),numel(nlambda));
PCOR=zeros(numel(alpha),numel(nlambda));
RSQ=zeros(numel(alpha),numel(nlambda));
for a=1:numel(alpha)
    alp=alpha(a);
    for nl=1:numel(nlambda)
        nlamb=nlambda(nl);
        opt=struct('alpha',alp,'nlambda',nlamb);
%         MSE=zeros(1,k);
        pred=zeros(size(InputMat,1),1);
        for kF=1:k
            idxTrain=training(c,kF);
            idxTest=test(c,kF);
            TrainInput=InputMat(idxTrain,:);
            TestInput=InputMat(idxTest,:);
            TrainObs=Observations(idxTrain,1);
            TestObs=Observations(idxTest,1);

            fit=cvglmnet(TrainInput,TrainObs,family,opt,[],10);
            pred(idxTest)=cvglmnetPredict(fit,TestInput);
            %MSE(1,kF)=immse(pred,TestObs); %MSE from default linear reg.
        end
%         avgMSE(a,nl)=mean(MSE);
%         stdMSE(a,nl)=std(MSE);
        [Corr,PCorr]=Correlation(pred,Observations);
        Rsq=Rsquared(pred,Observations);
        COR(a,nl)=Corr;
        PCOR(a,nl)=PCorr;
        RSQ(a,nl)=Rsq;
    end
end
clear a alp nl nlamb pred kF
CORRS=zeros(numel(alpha)*numel(nlambda),3);
PCORRS=zeros(numel(alpha)*numel(nlambda),3);
RSQS=zeros(numel(alpha)*numel(nlambda),3);
count=1;
for na=1:size(COR,1)
    alp=alpha(na);
    for nl=1:size(COR,2)
        nlamb=nlambda(nl);
        cor=COR(na,nl);
        pcor=PCOR(na,nl);
        rsq=RSQ(na,nl);
        CORRS(count,:)=[alp,nlamb,cor];
        PCORRS(count,:)=[alp,nlamb,pcor];
        RSQS(count,:)=[alp,nlamb,rsq];
        count=count+1;
    end
end
clear count nt np numT numP val
fileprefix=['./OptimizeGLMNet_' family];
fid = fopen([fileprefix '_Correlations.txt'],'wt');
fprintf(fid, '%s\t%s\t%s\n', 'Alpha','numLambda','Correlation');  % header
fclose(fid);
dlmwrite([fileprefix '_Correlations.txt'],CORRS,'-append','delimiter','\t','precision',10);

fid = fopen([fileprefix '_CorrelationPvalues.txt'],'wt');
fprintf(fid, '%s\t%s\t%s\n', 'Alpha','numLambda','Correlation_Pvalue');  % header
fclose(fid);
dlmwrite([fileprefix '_CorrelationPvalues.txt'],PCORRS,'-append','delimiter','\t','precision',10);

fid = fopen([fileprefix '_RSquared.txt'],'wt');
fprintf(fid, '%s\t%s\t%s\n', 'Alpha','numLambda','RSquared');  % header
fclose(fid);
dlmwrite([fileprefix '_RSquared.txt'],RSQS,'-append','delimiter','\t','precision',10);



end
