function [AUC, AUPR, REC1] = ROC(R,EFG,ExcludedRows);
    [m n] = size(R);
    no_excl = size(ExcludedRows,2);
    
    tot_hid = no_excl * n;
    ones_hid = sum(sum(R(ExcludedRows,:)));
    
    zeros_hid = tot_hid - ones_hid;
    
    ALL_EFG = zeros(m,n);
    ALL_EFG(:,:) = -1 * realmax;
    ALL_EFG(ExcludedRows,:) = EFG(ExcludedRows,:);
    
    [B Ind] = sort(ALL_EFG(:));
    [ii jj] = ind2sub([m n],Ind);
    DATAPT = flipud([ii jj]);
    DATAPT = DATAPT(1:tot_hid,:);
    DATAPT = flipud(DATAPT);
    
    TN = zeros(tot_hid + 1, 1);
    FN = zeros(tot_hid + 1, 1);
    
    TP = ones_hid * ones(tot_hid + 1, 1);
    FP = zeros_hid * ones(tot_hid + 1, 1);
    
    X_AUC(1) = 1; Y_AUC(1) = 1;
    X_AUPR(1) = 1; Y_AUPR(1) = 0;
    
    cntr = 2;
    
    for i=2:tot_hid
       x = DATAPT(i-1,1);
       y = DATAPT(i-1,2);
       if R(x,y) == 0
           TN(i) = TN(i-1) + 1;
           FP(i) = FP(i-1) - 1;
           TP(i) = TP(i-1);
           FN(i) = FN(i-1);
       else 
           TP(i) = TP(i-1) - 1;
           FN(i) = FN(i-1) + 1;           
           TN(i) = TN(i-1);
           FP(i) = FP(i-1);
       end
       
       xt = DATAPT(i,1);
       yt = DATAPT(i,2);
       if ALL_EFG(xt,yt) ~= ALL_EFG(x,y)

            X_AUC(cntr) = FP(i) / (FP(i) + TN(i));
            Y_AUC(cntr) = TP(i) / (TP(i) + FN(i));  
            
            X_AUPR(cntr) = TP(i) / (TP(i) + FN(i));
            Y_AUPR(cntr) = TP(i) / (TP(i) + FP(i));

            cntr = cntr + 1;
       end
    end
    
    i = tot_hid + 1;
    x = DATAPT(i-1,1);
    y = DATAPT(i-1,2);
    if R(x,y) == 0
        % X decreases, Y stays
        TN(i) = TN(i-1) + 1;
        FP(i) = FP(i-1) - 1;
        TP(i) = TP(i-1);
        FN(i) = FN(i-1);
    else 
       % X stays, Y decreases
       TP(i) = TP(i-1) - 1;
       FN(i) = FN(i-1) + 1;           
       TN(i) = TN(i-1);
       FP(i) = FP(i-1);
    end

    X_AUC(cntr) = FP(i) / (FP(i) + TN(i));
    Y_AUC(cntr) = TP(i) / (TP(i) + FN(i));  
    
    X_AUPR(cntr) = TP(i) / (TP(i) + FN(i));
    Y_AUPR(cntr) = 1;
    
    xauc = X_AUC(end:-1:1);
    yauc = Y_AUC(end:-1:1);
    
    xaupr = X_AUPR(end:-1:1);
    yaupr = Y_AUPR(end:-1:1);
    
    AUC = trapz(xauc, yauc);
    AUPR = trapz(xaupr, yaupr);
    % plot(xaupr,yaupr);
    % plot(xaupr,yaupr);
end

