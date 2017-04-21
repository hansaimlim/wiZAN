function glmstat()
load('C:\Users\hansaimlim\OneDrive - Cuny GradCenter\DrugMatrix Manual\AffyLiver\AffyLiverMeanResponse.mat')
corrs=zeros(1,10);
rsqs=zeros(1,10);
corrs_glm=zeros(1,10);
rsqs_glm=zeros(1,10);

for i=1:10
[pred,pred_glm]=GlmNetCoeffValidation;
c=corrcoef(AffyLiverMeanResponse(:,1),pred');
corrs(1,i)=c(1,2);
rsqs(1,i)=Rsquared(AffyLiverMeanResponse(:,1),pred');

cglm=corrcoef(AffyLiverMeanResponse(:,1),pred_glm');
corrs_glm(1,i)=cglm(1,2);
rsqs_glm(1,i)=Rsquared(AffyLiverMeanResponse(:,1),pred_glm');
end
disp(['coeff*gex: Mean correlation=' num2str(mean(corrs)) ' (std=' num2str(std(corrs)) ')' ...
    ' Mean R-squared=' num2str(mean(rsqs)) ' (std=' num2str(std(rsqs)) ')']);
disp(['glmnet: Mean correlation=' num2str(mean(corrs_glm)) ' (std=' num2str(std(corrs_glm)) ')' ...
    ' Mean R-squared=' num2str(mean(rsqs_glm)) ' (std=' num2str(std(rsqs_glm)) ')']);

end