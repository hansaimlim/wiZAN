function [ F ] = mulanImputeF( G,D_map,DM,alpha,beta,paras,F0 )
%MULANIMPUTE Summary of this function goes here
%   the function takes a multi-layered network as input and output the low
%   rank approximation for each layer
% \sum||W.(D(i,j)-Fi'Fj))||^2 +\sum \beta tr(Fi'(Ti-Ai)Fi)+ \sum beta||Fi||^2
% INPUT: G: adjacency matrix for each layers; D_map: dependency map matrix;
% DM: dependency matrice; beta:coefficient of ||Fi||; alpha:coefficients
% of user homopily;
% OUTPUT: set of low-rank approximations

%if no parameters were set
if nargin<6
    % weight,  rank, maxIte
    paras = [0.1,100, 100];
end
%if no initial low rank matrices are given, initialized them randomly
if nargin < 7
    rank = paras(2);
    for i = 1:length(G)
        m = size(G{i}.A,1);
        F0{i}.F =  (rand(m, rank))/sqrt(rank);
    end
end


[F] = updateF(G,D_map,DM,beta,alpha,F0,paras);

end

function [F] = updateF(G,D_map,DM,beta,alpha,F0,paras)
[I,J,V] = find(D_map);
weight = paras(1);
maxIter = paras(3);

threshold = 1e-8;

%initialized low rank approximations
for i = 1:length(G)
    F{i}.F = F0{i}.F;
end


maxChange = 1;
iter = 0;

%iteratively update each low-rank matrix
while maxChange > threshold && iter < maxIter
    %update each low-rank approximation matrix
    for i = 1:length(G)     
        i
        F0 = F{i}.F;
        F{i}.F = updateFi(i,D_map,weight,G,DM,F,beta,alpha);
        changes(i) = max(max(abs(F{i}.F-F0)));
    end
    maxChange  = max(changes)
    iter = iter + 1;
end
iter
end

function [Fi_new] = updateFi(curLayer,Dmap,weight,G,DM,F,beta,alpha)

Fi = F{curLayer}.F;
[m,n] = size(Fi);
%initialize A2 and B2
A2 = zeros(m,n);
B2 = zeros(m,n);
%find related layers
relLayers = find(Dmap(curLayer,:));
%get summation part of A2 and B2
for i = 1:length(relLayers)
    j = relLayers(i);
    index = Dmap(curLayer,j);
    Dij = DM{index}.D;
    if curLayer > relLayers(i)
        Dij = Dij';
    
    end
    A2 = A2 + Dij*F{j}.F;
    
    Rhat = get_UVT(Dij,Fi,F{j}.F);
    B2 = B2 + (1-weight)*Rhat*F{j}.F+weight*Fi*(F{j}.F'*F{j}.F);
end
%remaining part of A2
A2 = A2 + alpha*G{curLayer}.A*Fi;
%remaining parts of B2
sumA = sum(G{curLayer}.A,2);
DA = spdiags(sumA,0,m,m);
B2 = B2 + beta*Fi + alpha*DA*Fi + 1e-8;

Fi_new = Fi.*sqrt(A2./B2);

end


function [UVT] = get_UVT(R, U, V)

[m, n] = size(R);

[I, J] = find(R);
iSize = size(I, 1);
K = ones(iSize,1);
count = 0;
for j = 1:n
    id = find(R(:,j));
    len = length(id);
    K(count+1:count+len) = U(id, :) * V(j,:)';
    count = count + len;
end

UVT = sparse(I, J, K, m, n);

end