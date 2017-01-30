function [TPR time] = TprBench(R,M,N,T,J,lR,lM,lN,iter,rnk,WP,M_cut,N_cut,QUICK)
    tic
    
    % to perform TPR analysis on drugs, 
    % we need to switch things around
    
    % if CVS == 2
    	R = R';
    	M1 = M';
    	N1 = N';
    	M = N1;
    	N = M1;
    	lM1 = lM;
    	lM = lN;
    	lN = lM1;
    	clear M1
    	clear N1
    % end     
    
    TANIM = 0.7; FGSCORE = 0.95; HWGHT = 7; LWGHT = 2; WGHT = 7;

    [m n] = size(R)
    
    [DM nM]= GetDiag(M,J);
    [DN nN]= GetDiag(N,J);

    DMM = DM-nM;
    DNN = DN-nN;        
    
    [p q] = size(T);

    TEST_MTX = R;

    for k=1:p
        i = T(k,1);
        j = T(k,2);
        TEST_MTX(i,j) = 0;
    end
    
    [ExcludedRows unimportant] = find(sum(TEST_MTX,2)==0); 
    [unimportant ExcludedColumns] = find(sum(TEST_MTX,1)==0); 
                
    W = max(1, 6 * TEST_MTX);
    Q = zeros(m,n);            

    [F G] = WeightImputeLogFactorization(TEST_MTX,DMM,DNN,W,Q,lR,lM,lN,iter,rnk);

    [nF nG HI_IND] = WeightedProfile(F, G, nM, nN, ExcludedRows, ExcludedColumns, J + 2, TANIM, WP, M_cut, N_cut);

    EXC = GetP(nF*nG'); 
    if QUICK == 0
        
%         LOW = EXC < FGSCORE;
%         TOO_LOW = EXC < 1-FGSCORE;
%         HIGH = EXC >= FGSCORE;
%         EXC(HIGH) = 1;
%         EXC(LOW) = 0;
%         W(HIGH) = HWGHT;
%         W(TOO_LOW) = LWGHT;
%         EXC = max(EXC,TEST_MTX);
% 
%         [F G] =  WeightImputeLogFactorization(EXC,DMM,DNN,W,Q,lR,lM,lN,iter,rnk);  
%     
%         [nF nG] = WeightedProfile(F, G, nM, nN, ExcludedRows, ExcludedColumns, J + 2, TANIM, WP, M_cut, N_cut);
%         EXC = GetP(nF*nG'); 

        LOW = EXC < FGSCORE;
        HIGH = EXC >= FGSCORE;
        EXC(HIGH) = 1;
        EXC(LOW) = 0;
        W(HI_IND > 0, :) = WGHT;
        EXC = max(EXC,TEST_MTX);

        [F G] =  WeightImputeLogFactorization(EXC,DMM,DNN,W,Q,lR,lM,lN,iter,rnk);  

        [nF nG HI_IND] = WeightedProfile(F, G, nM, nN, ExcludedRows, ExcludedColumns, J + 2, TANIM, WP, M_cut, N_cut);
        EXC = GetP(nF*nG'); 
    end
    
    TPR = ComputeTPR(T, EXC);

    time = toc;
end

function [ TPR ] = ComputeTPR(T, MTX)
    [m n] = size(MTX);
    t1p = ceil(n / 100);
    cntr = 0;
    [numpairs ignore] = size(T);
    for k=1:numpairs
        i = T(k,1);
        j = T(k,2);
        SORTED_ROW = sort(MTX(i,:),'descend');
        if MTX(i,j) >= SORTED_ROW(t1p) 
            cntr = cntr+1;
        end
    end
    TPR = cntr / numpairs;
end

