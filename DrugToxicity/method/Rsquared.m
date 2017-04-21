function Rsq=Rsquared(pred,Observations)
    sse=sum((pred-Observations).^2);
    sst=sum( (Observations-mean(Observations)).^2 );
    Rsq=1-(sse/sst);
end