function [tF tG tempS] = WeightImputeLogFactorization(R,DMM,DNN,W,Q,lR,lM,lN,iter,rnk)
    % Note: learn_rate is the initial learning rate in AdaGrad,
    % which can be set to 0.5 for many problems.
    
    rng shuffle;
    
    [m, n] = size(R);
    
    learn_rate = 0.5;

    tF = normrnd(0,1/sqrt(rnk),[m,rnk]);
    tG = normrnd(0,1/sqrt(rnk),[n,rnk]);
    
    tGt = (tG)';
    
    tempS = 0;
    
    Im = eye(m);
    In = eye(n);

    df_sum = zeros(m,rnk);
    dg_sum = zeros(n,rnk);
       
    itr = 0;
    
    Wt = W';

    RpQ = R+Q;
    RpQt = (RpQ)';
      
    while itr < iter
                        
        itr = itr + 1;

        P = GetP(tF*tGt);
        df = (W .* (P - (RpQ))) * tG + lR.*tF + lM.*(DMM*tF);
        df_sum = df_sum + (df .* df);

        tF = tF - df .* (learn_rate ./ sqrt(df_sum));

        P = GetP(tF*tGt);
        Pt=P';
        
        dg = (Wt .* (Pt - RpQt)) * tF + lR.*tG + lN.*(DNN*tG);
        dg_sum = dg_sum + (dg .* dg);
        
        tG = tG - dg .* (learn_rate ./ sqrt(dg_sum));
        tGt = (tG)';
    end
    % fprintf('\n');
    % tFt = (tF)';        
    % FGt = tF * tGt;
    % B = W .* (log(1+exp(FGt)) - (RpQ.*FGt));
    % X = sum(B(:));
    % Y = 0.5 * trace(tFt*(lR.*Im + lM.*DMM)*tF);
    % Z = 0.5 * trace(tGt*(lR.*In + lN.*DNN)*tG); 
    % tempS = X + Y + Z;
end
