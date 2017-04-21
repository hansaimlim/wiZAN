function SVMBenchmark(InputMat,Observations,kFold,KernelFunction)
maxNumCompThreads(3); %determine the maximum number of cores to use
k=kFold;

%Compare LeastSquares with SVM
[lrmdl,FitInfo,HyperParamOptRes]=fitrlinear(InputMat,Observations,...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'));
%find optimized hyperparameters
svmgau=fitrsvm(InputMat,Observations,'Standardize',true,...
    'KernelFunction',KernelFunction,'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'));


bestBC=mean(svmgau.BoxConstraints);
bestKS=svmgau.KernelParameters.Scale;
bestEps=svmgau.Epsilon;
disp(['Best parameters: BoxConstraints=' num2str(bestBC) ', KernelScale=' ...
    num2str(bestKS) ', Epsilon=' num2str(bestEps)]);
c=cvpartition(size(InputMat,1),'kFold',k);
pred=zeros(size(InputMat,1),1);
pred_lsq=zeros(size(InputMat,1),1);
for kF=1:k
    idxTrain=training(c,kF);
    idxTest=test(c,kF);
    TrainInput=InputMat(idxTrain,:);
    TestInput=InputMat(idxTest,:);
    TrainObs=Observations(idxTrain,1);
    TestObs=Observations(idxTest,1);

    SVMGAU=fitrsvm(TrainInput,TrainObs,'KernelFunction',KernelFunction,...
    'KernelScale',bestKS,'Epsilon',bestEps,'BoxConstraint',bestBC,...
    'Standardize',true);
    pred(idxTest)=predict(SVMGAU,TestInput);
    LSQ=fitrlinear(TrainInput,TrainObs,'Learner','leastsquares');
    pred_lsq(idxTest)=predict(LSQ,TestInput);
    %MSE(1,kF)=immse(pred,TestObs); %MSE from default linear reg.
end

rsqSVMgau=Rsquared(pred,Observations);
[corSVMgau,pval]=Correlation(pred,Observations);
disp(['SVM ' KernelFunction ' Cross validation. R-squared=' num2str(rsqSVMgau) ...
    ', Correlation=' num2str(corSVMgau)]);

rsqLSQ=Rsquared(pred_lsq,Observations);
[corLSQ,pval_lsq]=Correlation(pred_lsq,Observations);
disp(['Least Squares Cross validation. R-squared=' num2str(rsqLSQ) ...
    ', Correlation=' num2str(corLSQ)]);



end
