function [nF, nG, HI_IND] = WeightedProfile(F, G, M, N, ExcludedRows, ExcludedColumns, J, TANIM, WP, M_cut, N_cut)

    m = size(M,1);
    n = size(N,1);
        
    HI_IND = zeros(1,m);
    
    nF = F;
    for i = 1:m
       if sum(any(i==ExcludedRows))>0
            Row = M(i,:);       
            Row(ExcludedRows) = -1 * realmax;
            summation = 0;
            nF(i,:) = WP * F(i,:);
            for j=1:J
                [Mx, Ix] = max(Row);
                if Mx > TANIM
                   HI_IND(i) = 1;
                end
                if Mx > M_cut % M_cut = 0.3 for ZINC and -1 otherwise
                   summation = summation + Mx;
                   nF(i,:) = nF(i,:) + Mx * F(Ix,:);
                end
                Row(Ix) = -1 * realmax;
            end
            nF(i,:) = nF(i,:) / (WP + summation);
       end
       
    end
    
    nG = G;
    for i = 1:n
       if sum(any(i==ExcludedColumns))>0
            Col = N(i,:);
            Col(ExcludedColumns) = -1 * realmax;
            summation = 0;
            nG(i,:) = WP * G(i,:);
            for j=1:J
                [Nx Ix] = max(Col);
                if Nx > N_cut % N_cut = 0.5 for ZINC and -1 otherwise
                  summation = summation + Nx;
                  nG(i,:) = nG(i,:) + Nx * G(Ix,:);
                end
                Col(Ix) = -1 * realmax;
            end
            if N_cut > 0
                nG(i,:) = nG(i,:) / (WP + summation); 
            end
       end
    end 
end
