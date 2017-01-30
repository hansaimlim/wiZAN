function P=COSINE_Predict(R,M,N,rnk,lR,lM,lN,WP,imp,FacOpt)
% Using user-defined parameters to get prediction matrix
% Returns (m by n) matrix P, containing predicted scores for each pair

% R: knwon association matrix (m by n)
% M: similarity score matrix (m by m)
% N: similarity score matrix (n by n)
% rnk: Rank parameter (integer, maximum_rank=min(m,n) )
% lR: lambda_R (scalar in [0,1])
% lM: lambda_M (scalar in [0,1])
% lN: lambda_N (scalar in [0,1])
% WP: weight (scalar in [0,1])
% imp: imputation value (scalar in [0,1])
% FacOpt: 'Lin' for Linear Factorization; anything else for Logistic Factorization
if lR>1
    msg='Parameter lR should be between 0 and 1.';
    error(msg);
elseif lM>1
    msg='Parameter lM should be between 0 and 1.';
    error(msg);
elseif lN>1
    msg='Parameter lN should be between 0 and 1.';
    error(msg);
elseif WP>1
    msg='Parameter WP should be between 0 and 1.';
    error(msg);
elseif imp>1
    msg='Parameter imp should be between 0 and 1.';
    error(msg);
end
[m,n]=size(R);
if size(M,1)~=m
    msg='Matrix size does not match.';
    error(msg);
elseif size(N,1)~=n
    msg='Matrix size does not match.';
    error(msg);
end
if rnk > min(m,n)
    msg=['Rank parameter is too large. Maximum rank is ' num2str(min(m,n))...
        ' for your data.'];
    error(msg);
end
if strcmp(FacOpt, 'Lin')
    disp('Linear Factorization selected.');
else
    disp('Logistic Factorization selected.'); 
end
    W = max(1, 6 * R);
    Q = ~R.*imp;    
    iter=100; %default iteration number
    TANIM = 0.7; FGSCORE = 0.95; WGHT = 7; J = 2;
    M_cut = -1;
    N_cut = -1;
    [DM, nM]= GetDiag(M,J);
    [DN, nN]= GetDiag(N,J);
    DMM = DM-nM;
    DNN = DN-nN;
    [~,ColdStartCols]=find(sum(R,1)==0);
    [ColdStartRows,~]=find(sum(R,2)==0);
    if strcmp(FacOpt, 'Lin')
        [F, G] = WeightImputeLinFactorization(R,M,N,W,Q,lR,lM,lN,iter,rnk);
        [nF, nG, HI_IND] = WeightedProfile(F, G, nM, nN, ColdStartRows, ColdStartCols, J+2, TANIM, WP, M_cut, N_cut);
        EXC = nF*nG';
    else
        [F, G] = WeightImputeLogFactorization(R,DMM,DNN,W,Q,lR,lM,lN,iter,rnk);
        [nF, nG, HI_IND] = WeightedProfile(F, G, nM, nN, ColdStartRows, ColdStartCols, J+2, TANIM, WP, M_cut, N_cut);
        EXC = GetP(nF*nG'); 
    end

    LOW = EXC < FGSCORE;
    HIGH = EXC >= FGSCORE;
    EXC(HIGH) = 1;
    EXC(LOW) = 0;
    W(HI_IND > 0, :) = WGHT;
    EXC = max(EXC,R);

    if strcmp(FacOpt, 'Lin')
        [F, G] =  WeightImputeLinFactorization(EXC,M,N,W,Q,lR,lM,lN,iter,rnk);    
        [nF, nG, ~] = WeightedProfile(F, G, nM, nN, ColdStartRows, ColdStartCols, J+2, TANIM, WP, M_cut, N_cut);
    else
        [F, G] =  WeightImputeLogFactorization(EXC,DMM,DNN,W,Q,lR,lM,lN,iter,rnk);  
        [nF, nG, ~] = WeightedProfile(F, G, nM, nN, ColdStartRows, ColdStartCols, J+2, TANIM, WP, M_cut, N_cut);
    end
    
    
    if strcmp(FacOpt, 'Lin')
        P=nF*nG';
    else
        P=GetP(nF*nG'); 
    end

end