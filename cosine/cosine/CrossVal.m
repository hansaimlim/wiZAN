function [AUC_AVG AUC_CI AUPR_AVG AUPR_CI time] = CrossVal(R,M,N,J,lR,lM,lN,iter,rnk,folds,rounds,CV2,LINFAC,WP,M_cut,N_cut)
    tic
    
    TANIM = 0.7; FGSCORE = 0.95; WGHT = 7;
    
    if CV2 == 1
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
    end

    [m n] = size(R);
    
    [DM nM]= GetDiag(M,J);
    [DN nN]= GetDiag(N,J);

    DMM = DM-nM;
    DNN = DN-nN;        
    
    AUC = zeros(rounds,1);
    AUPR = zeros(rounds,1);

    if mod(m,folds) == 0
        binsize = m/folds - 1;  
    else
        binsize = floor(m/folds);
    end
    
    for j=1:rounds
        RowPerm = randperm(m);
        
        auc_sum = 0;
        apr_sum = 0;
        
        for t = 1:folds
            
            rng shuffle;
            
            ExcludedRows = sort(RowPerm((t-1)*binsize+1:t*binsize));
            if t <= m-binsize*folds
                ExcludedRows = sort([ExcludedRows RowPerm(t+binsize*folds)]);
            end
            
            TEST_MTX = R;

            TEST_MTX(ExcludedRows,:) = 0;
            [unimportant ExcludedColumns] = find(sum(TEST_MTX,1)==0); 
            
            W = max(1, 6 * TEST_MTX);
            Q = zeros(m,n);            
                        
            if LINFAC == 1
                [F G] = WeightImputeLinFactorization(TEST_MTX,M,N,W,Q,lR,lM,lN,iter,rnk);
                [nF nG HI_IND] = WeightedProfile(F, G, nM, nN, ExcludedRows, ExcludedColumns, J+2, TANIM, WP, M_cut, N_cut);
                EXC = nF*nG';
            else
                [F G] = WeightImputeLogFactorization(TEST_MTX,DMM,DNN,W,Q,lR,lM,lN,iter,rnk);
                [nF nG HI_IND] = WeightedProfile(F, G, nM, nN, ExcludedRows, ExcludedColumns, J+2, TANIM, WP, M_cut, N_cut);
                EXC = GetP(nF*nG'); 
            end
            
%         LOW = EXC < FGSCORE;
%         TOO_LOW = EXC < 1-FGSCORE;
%         HIGH = EXC >= FGSCORE;
%         EXC(HIGH) = 1;
%         EXC(LOW) = 0;
%         W(HIGH) = HWGHT;
%         W(TOO_LOW) = LWGHT;
%         EXC = max(EXC,TEST_MTX);

            LOW = EXC < FGSCORE;
            HIGH = EXC >= FGSCORE;
            EXC(HIGH) = 1;
            EXC(LOW) = 0;
            W(HI_IND > 0, :) = WGHT;
            EXC = max(EXC,TEST_MTX);
                        
            if LINFAC == 1
                [F G] =  WeightImputeLinFactorization(EXC,M,N,W,Q,lR,lM,lN,iter,rnk);    
                [nF nG HI_IND] = WeightedProfile(F, G, nM, nN, ExcludedRows, ExcludedColumns, J+2, TANIM, WP, M_cut, N_cut);
                [ac apr] = ROC(R,nF*nG',ExcludedRows);
            else
                [F G] =  WeightImputeLogFactorization(EXC,DMM,DNN,W,Q,lR,lM,lN,iter,rnk);  
                [nF nG HI_IND] = WeightedProfile(F, G, nM, nN, ExcludedRows, ExcludedColumns, J+2, TANIM, WP, M_cut, N_cut);
                [ac apr] = ROC(R,GetP(nF*nG'),ExcludedRows);
            end
            
            auc_sum = auc_sum + ac;
            apr_sum = apr_sum + apr;
                       
        end
        
        AUC(j) = auc_sum / folds;
        AUPR(j) = apr_sum / folds;
    end
    
    [AUC_AVG AUC_CI] = Confidence(AUC,0.95);
    [AUPR_AVG AUPR_CI] = Confidence(AUPR,0.95);
    time = toc;
end

