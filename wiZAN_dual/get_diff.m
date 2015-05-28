function [MAP, MPR, HLU, AUC] = get_diff(test, U, V, para)
%{
real = test(:,3);
t_size = size(test, 1);
predict = zeros(t_size,1);

for i = 1:t_size
    x = test(i,1);
    y = test(i,2);
    predict(i,1) = U(x,:) * V(y,:)';
end


topN = (1000:100:2000);
[AUC, avgF, avgP, avgR] = computeAUC(real, predict, topN);
%}

m = size(U, 1);
n = size(V, 1);

% para: lambda, r, T, rank, maxIte, ite_of_bisection method, topN, gamma, lambda

test = sparse(test(:,1), test(:,2), test(:,3), m, n);

MAP = zeros(m,1);
HLU_1 = zeros(m,1);
HLU_2 = zeros(m,1);
PR_1 = zeros(m,1);
PR_2 = zeros(m,1);
AUC = zeros(m,1);

for i = 1:m
    real = test(i,:)';
    predict = (U(i,:) * V')';
    
    MAP(i,1) = computeAP(real, predict);
    [HLU_1(i,1), HLU_2(i,1)] = computeHLU(real, predict);
    [PR_1(i,1), PR_2(i,1)] = computePR(real, predict);
    AUC(i,1) = computeAUC(real, predict);
end

[I] = find(MAP>0);
MAP = MAP(I);
MAP = mean(MAP);
HLU = 100 * sum(HLU_1) / sum(HLU_2);
MPR = sum(PR_1) / sum(PR_2);
[I] = find(AUC>0);
AUC = AUC(I);
AUC = mean(AUC);

end








