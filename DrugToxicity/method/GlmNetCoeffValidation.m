function [preds,preds_glm]=GlmNetCoeffValidation()
%Uses coefficients generated using GlmNet 
%to predict values in observations by regression
%Uses optimized parameters of alpha and nlambda for AffyLiver dataset
maxNumCompThreads(3);
k=10; %kFold cross validation
addpath('C:\Users\hansaimlim\OneDrive - Cuny GradCenter\DrugMatrix Manual\algo\LASSO\glmnet_matlab');
load('C:\Users\hansaimlim\OneDrive - Cuny GradCenter\DrugMatrix Manual\AffyLiver\AffyLiver.mat');
load('C:\Users\hansaimlim\OneDrive - Cuny GradCenter\DrugMatrix Manual\AffyLiver\AffyLiverMeanResponse.mat');
InputMat=AffyLiver;
Observations=AffyLiverMeanResponse;

opt=struct('alpha',0.8,'nlambda',829); %optimized params for AffyLiver
c=cvpartition(size(InputMat,1),'kFold',k);
preds=zeros(1,size(InputMat,1));
cors=zeros(1,k);
rsqs=zeros(1,k);

preds_glm=zeros(1,size(InputMat,1));
cors_glm=zeros(1,k);
rsqs_glm=zeros(1,k);

for kF=1:k
    idxTrain=training(c,kF);
    idxTest=test(c,kF);
    TrainInput=InputMat(idxTrain,:);
    TestInput=InputMat(idxTest,:);
    TrainObs=Observations(idxTrain,1);
    TestObs=Observations(idxTest,1);

    glmgau=glmnet(TrainInput,TrainObs,'gaussian',opt);
    glmgaubeta=glmgau.beta;
    glmgaucoeff=glmgaubeta(:,end);
    pred=TestInput*glmgaucoeff;
    preds(idxTest)=pred;
    [cor,pcor]=Correlation(pred,TestObs);
    rsq=Rsquared(pred,TestObs);
    cors(1,kF)=cor;
    rsqs(1,kF)=rsq;
    
    fit=cvglmnet(TrainInput,TrainObs,'gaussian',opt,[],10);
    predglm=cvglmnetPredict(fit,TestInput);
    preds_glm(idxTest)=predglm;
    [cor_glm,pcor_glm]=Correlation(predglm,TestObs);
    rsq_glm=Rsquared(predglm,TestObs);
    cors_glm(1,kF)=cor_glm;
    rsqs_glm(1,kF)=rsq_glm;
    %MSE(1,kF)=immse(pred,TestObs); %MSE from default linear reg.
end
cors(isnan(cors))=0;
rsqs(isnan(rsqs))=0;
cors_glm(isnan(cors_glm))=0;
rsqs_glm(isnan(rsqs_glm))=0;
disp(['(Coefficient*g-ex):10Fold Mean Correlation=' num2str(mean(cors)) ' std(' num2str(std(cors)) ')' ...
    '; Mean R-squared=' num2str(mean(rsqs)) ' std(' num2str(std(rsqs)) ')']);
disp(['GlmNetGaussian:10Fold Mean Correlation=' num2str(mean(cors_glm)) ' std(' num2str(std(cors_glm)) ')' ...
    '; Mean R-squared=' num2str(mean(rsqs_glm)) ' std(' num2str(std(rsqs_glm)) ')']);


end