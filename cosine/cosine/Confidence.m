function [mn ci] = Confidence(x,conf)
% use conf interval for NRLMF benchmark
    pd = fitdist(x,'Normal');
    ci = paramci(pd);
    mn = mean(pd);
    ci = ci(2) - mn;

% use std for Laarhoven benchmark
% mn = mean(x);
% ci = std(x);
end

