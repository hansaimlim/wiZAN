function [] = Run()
    % add path to files here
    % addpath('');

    TARG = 'Nr' % or 'Gpcr' or 'Ion' or 'Enz'
    LINFAC = 0 % results in logistic factorization; otherwise linear
    CV2 = 1; % is for "new drugs"; otherwise new proteins

    rng shuffle;

    folds = 10; rounds = 5; J = 5; 

    fn = strcat(TARG,'_R08');
    load ([fn '.dat']);
    R = spconvert(eval(fn));
    fn = strcat(TARG,'_M08');
    load ([fn '.dat']);
    M = spconvert(eval(fn));
    fn = strcat(TARG,'_N08');
    load ([fn '.dat']);
    N = spconvert(eval(fn));

    [m,n] = size(R);

    maxiter = floor(m * n * 0.001591425 + 30.23);
    iter = min(maxiter,100);
    rnk = min(m,n);

    fprintf('%s FOLDS:%d ROUNDS:%d\n',TARG,folds,rounds);
    
    WP = 0.0;
    M_cut = -1;
    N_cut = -1;
    
    % need to train to get optimal parameters for different classes
    lR=0.01; lM=1.0; lN=0.1;
    
    [AUC_AVG AUC_CI AUPR_AVG AUPR_CI time] = CrossVal(R,M,N,J,lR,lM,lN,iter,rnk,folds,rounds,CV2,LINFAC,WP,M_cut,N_cut);
    fprintf('AUC:%f AUPR:%f lR:%f lM:%f lN:%f time:%f maxiter:%d\n',AUC_AVG,AUPR_AVG,lR,lM,lN,time,iter);
end
