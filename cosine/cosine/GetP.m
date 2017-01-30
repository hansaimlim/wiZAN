function [ P ] = GetP(X)
    P = exp(X);
    P = P ./ (1 + P);
end