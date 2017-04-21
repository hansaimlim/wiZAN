function [cor,pval]=Correlation(pred,Observations)
    [corMat,pvalMat]=corrcoef(pred,Observations);
    cor=corMat(1,2);
    pval=pvalMat(1,2);
end