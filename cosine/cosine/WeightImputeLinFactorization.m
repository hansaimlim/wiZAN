function [tF, tG, tempS] = WeightImputeLinFactorization(R,M,N,W,Q,lR,lM,lN,iter,rnk)

    rng shuffle;
    
    [i_row, j_col] = find(R);
    [m, n] = size(R);

    DM = diag(sum(M,2));
    DN = diag(sum(N,2));
    
    tF = rand(m, rnk);
    tG = rand(n, rnk);
    
    frobF_new = realmax;
    frobG_new = realmax;

    Rt = R';
    RpQ = R+Q;
    RpQt = (RpQ)';
    
    MINDIFF = 0.00001;
    diff = realmax;
    itr = 0;
 
    while diff > MINDIFF && itr < iter
        
        itr = itr + 1;
        
        frobF = frobF_new;
        frobG = frobG_new;

        R1 = sparse(m,n);
        for x = 1:numel(i_row)
            u = i_row(x);
            i = j_col(x);
            R1(u,i) = tF(u,:)*tG(i,:)';
        end
   
        A1 = (W.*W.*RpQ)*tG + lM*M*tF;
        B1 = (W.*W.*(tF*tG'))*tG + lR*tF + lM*DM*tF;

        A2 = (W'.*W'.*RpQt)*tF + lN*N*tG;
        B2 = (W'.*W'.*(tG*tF'))*tF + lR*tG + lN*DN*tG;

        Fincrement = sqrt(A1./B1);
        tF = Fincrement.*tF;
    
        Gincrement = sqrt(A2./B2);
        tG = Gincrement.*tG;
    
        frobF_new = norm(tF,'fro');
        frobG_new = norm(tG,'fro');
    
        diff = max(abs(frobF-frobF_new),abs(frobG-frobG_new));
        
%         tFG = norm(W.*(RpQ-tF*tG'),'fro').^2;
%         tLf = lM*trace(tF'*(DM-M)*tF);
%         tLg = lN*trace(tG'*(DN-N)*tG);
%         tR = lR*(norm(tF,'fro').^2 + norm(tG,'fro').^2);         
%         tempS = tFG + tLf + tLg + tR;
         
%         if tempS < S
%             S = tempS;
%             F = tF;
%             G = tG;
%         end   
    end
end
