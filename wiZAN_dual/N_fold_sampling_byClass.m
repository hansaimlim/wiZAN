function N_fold_sampling_byClass(chem_prot, chem_chem, prot_prot)
%returns training and test sets from input R matrix
%each training+test equals to the R (e.g. train1 + test1 = train2 + test2 = R)
%N=10;	%N fold
%rline='some_csv_file.csv';
R=chem_prot;
numedge=sum(R(:)>0);
nchem=size(R,1);
nprot=size(R,2);

chem_chem=chem_chem-speye(size(chem_chem)); %remove self-identity
prot_prot=prot_prot-speye(size(prot_prot)); %remove self-identity

%classify functions
%NT=sum(chem_prot(chemIdx,:)); %get number of targets by chemical index
%NL=sum(chem_prot(:,protIdx)); %get number of ligands by protein index
%maxTC=max(chem_chem(chemIdx,:));
%classification end

N3upL21more_rc=[];
N3upL16to20_rc=[];
N3upL11to15_rc=[];
N3upL6to10_rc=[];
N3upL1to5_rc=[];

N2L21more_rc=[];
N2L16to20_rc=[];
N2L11to15_rc=[];
N2L6to10_rc=[];
N2L1to5_rc=[];

N1L21more_rc=[];
N1L16to20_rc=[];
N1L11to15_rc=[];
N1L6to10_rc=[];
N1L1to5_rc=[];

N3upTc9to10_rc=[];
N3upTc8to9_rc=[];
N3upTc7to8_rc=[];
N3upTc6to7_rc=[];
N3upTc5to6_rc=[];

N2Tc9to10_rc=[];
N2Tc8to9_rc=[];
N2Tc7to8_rc=[];
N2Tc6to7_rc=[];
N2Tc5to6_rc=[];

N1Tc9to10_rc=[];
N1Tc8to9_rc=[];
N1Tc7to8_rc=[];
N1Tc6to7_rc=[];
N1Tc5to6_rc=[];
[chemIdxs,protIdxs]=find(chem_prot);
for i=1:numel(chemIdxs)
    chemIdx=chemIdxs(i);
    protIdx=protIdxs(i);
    NT=sum(chem_prot(chemIdx,:)); %get number of targets by chemical index
    NL=sum(chem_prot(:,protIdx)); %get number of ligands by protein index
    maxTC=max(chem_chem(chemIdx,:)); %get maximum chemical-chemical similarity
    if NT>=3
        %N3up
        if NL>=21
            %NL21more
            N3upL21more_rc=[N3upL21more_rc;chemIdx,protIdx];
        elseif NL>=16
            %NL16to20
            N3upL16to20_rc=[N3upL16to20_rc;chemIdx,protIdx];
        elseif NL>=11
            %NL11to15
            N3upL11to15_rc=[N3upL11to15_rc;chemIdx,protIdx];
        elseif NL>=6
            %NL6to10
            N3upL6to10_rc=[N3upL6to10_rc;chemIdx,protIdx];
        elseif NL>=1
            %NL1to5
            N3upL1to5_rc=[N3upL1to5_rc;chemIdx,protIdx];
        else
            %protein has no known ligands
        end
       
        if maxTC>=0.9
            %Tc0.9to1
            N3upTc9to10_rc=[N3upTc9to10_rc;chemIdx,protIdx];
        elseif maxTC>=0.8
            %Tc0.8to0.9
            N3upTc8to9_rc=[N3upTc8to9_rc;chemIdx,protIdx];
        elseif maxTC>=0.7
            %Tc0.7to0.8
            N3upTc7to8_rc=[N3upTc7to8_rc;chemIdx,protIdx];
        elseif maxTC>=0.6
            %Tc0.6to0.7
            N3upTc6to7_rc=[N3upTc6to7_rc;chemIdx,protIdx];
        elseif maxTC>=0.5
            %Tc0.5to0.6
            N3upTc5to6_rc=[N3upTc5to6_rc;chemIdx,protIdx];
        else
            %no similar chemicals
        end
    elseif NT==2
        %N2
        if NL>=21
            %NL21more
            N2L21more_rc=[N2L21more_rc;chemIdx,protIdx];
        elseif NL>=16
            %NL16to20
            N2L16to20_rc=[N2L16to20_rc;chemIdx,protIdx];
        elseif NL>=11
            %NL11to15
            N2L11to15_rc=[N2L11to15_rc;chemIdx,protIdx];
        elseif NL>=6
            %NL6to10
            N2L6to10_rc=[N2L6to10_rc;chemIdx,protIdx];
        elseif NL>=1
            %NL1to5
            N2L1to5_rc=[N2L1to5_rc;chemIdx,protIdx];
        else
            %protein has no known ligands
        end
       
        if maxTC>=0.9
            %Tc0.9to1
            N2Tc9to10_rc=[N2Tc9to10_rc;chemIdx,protIdx];
        elseif maxTC>=0.8
            %Tc0.8to0.9
            N2Tc8to9_rc=[N2Tc8to9_rc;chemIdx,protIdx];
        elseif maxTC>=0.7
            %Tc0.7to0.8
            N2Tc7to8_rc=[N2Tc7to8_rc;chemIdx,protIdx];
        elseif maxTC>=0.6
            %Tc0.6to0.7
            N2Tc6to7_rc=[N2Tc6to7_rc;chemIdx,protIdx];
        elseif maxTC>=0.5
            %Tc0.5to0.6
            N2Tc5to6_rc=[N2Tc5to6_rc;chemIdx,protIdx];
        else
            %no similar chemicals
        end
    elseif NT==1
        %N1
        if NL>=21
            %NL21more
            N1L21more_rc=[N1L21more_rc;chemIdx,protIdx];
        elseif NL>=16
            %NL16to20
            N1L16to20_rc=[N1L16to20_rc;chemIdx,protIdx];
        elseif NL>=11
            %NL11to15
            N1L11to15_rc=[N1L11to15_rc;chemIdx,protIdx];
        elseif NL>=6
            %NL6to10
            N1L6to10_rc=[N1L6to10_rc;chemIdx,protIdx];
        elseif NL>=1
            %NL1to5
            N1L1to5_rc=[N1L1to5_rc;chemIdx,protIdx];
        else
            %protein has no known ligands
        end
       
        if maxTC>=0.9
            %Tc0.9to1
            N1Tc9to10_rc=[N1Tc9to10_rc;chemIdx,protIdx];
        elseif maxTC>=0.8
            %Tc0.8to0.9
            N1Tc8to9_rc=[N1Tc8to9_rc;chemIdx,protIdx];
        elseif maxTC>=0.7
            %Tc0.7to0.8
            N1Tc7to8_rc=[N1Tc7to8_rc;chemIdx,protIdx];
        elseif maxTC>=0.6
            %Tc0.6to0.7
            N1Tc6to7_rc=[N1Tc6to7_rc;chemIdx,protIdx];
        elseif maxTC>=0.5
            %Tc0.5to0.6
            N1Tc5to6_rc=[N1Tc5to6_rc;chemIdx,protIdx];
        else
            %no similar chemicals
        end
    else
        %0 known target. Cold-start chemicals
    end
    
    
end

%N3L21more=sparse(N3upL21more_rc(:,1),N3upL21more_rc(:,2),1,nchem,nprot);
%N3L1620=sparse(N3upL16to20_rc(:,1),N3upL16to20_rc(:,2),1,nchem,nprot);
%N3L1115=sparse(N3upL11to15_rc(:,1),N3upL11to15_rc(:,2),1,nchem,nprot);
%N3L610=sparse(N3upL6to10_rc(:,1),N3upL6to10_rc(:,2),1,nchem,nprot);
%N3L15=sparse(N3upL1to5_rc(:,1),N3upL1to5_rc(:,2),1,nchem,nprot);
%N_fold_sampling(R, N3L21more, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTNL/N3/L21more/');
%N_fold_sampling(R, N3L1620, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTNL/N3/L16to20/');
%N_fold_sampling(R, N3L1115, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTNL/N3/L11to15/');
%N_fold_sampling(R, N3L610, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTNL/N3/L6to10/');
%N_fold_sampling(R, N3L15, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTNL/N3/L1to5/');

%N2L21more=sparse(N2L21more_rc(:,1),N2L21more_rc(:,2),1,nchem,nprot);
%N2L1620=sparse(N2L16to20_rc(:,1),N2L16to20_rc(:,2),1,nchem,nprot);
%N2L1115=sparse(N2L11to15_rc(:,1),N2L11to15_rc(:,2),1,nchem,nprot);
%N2L610=sparse(N2L6to10_rc(:,1),N2L6to10_rc(:,2),1,nchem,nprot);
%N2L15=sparse(N2L1to5_rc(:,1),N2L1to5_rc(:,2),1,nchem,nprot);
%N_fold_sampling(R, N2L21more, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTNL/N2/L21more/');
%N_fold_sampling(R, N2L1620, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTNL/N2/L16to20/');
%N_fold_sampling(R, N2L1115, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTNL/N2/L11to15/');
%N_fold_sampling(R, N2L610, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTNL/N2/L6to10/');
%N_fold_sampling(R, N2L15, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTNL/N2/L1to5/');

%N1L21more=sparse(N1L21more_rc(:,1),N1L21more_rc(:,2),1,nchem,nprot);
%N1L1620=sparse(N1L16to20_rc(:,1),N1L16to20_rc(:,2),1,nchem,nprot);
%N1L1115=sparse(N1L11to15_rc(:,1),N1L11to15_rc(:,2),1,nchem,nprot);
%N1L610=sparse(N1L6to10_rc(:,1),N1L6to10_rc(:,2),1,nchem,nprot);
%N1L15=sparse(N1L1to5_rc(:,1),N1L1to5_rc(:,2),1,nchem,nprot);
%N_fold_sampling(R, N1L21more, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTNL/N1/L21more/');
%N_fold_sampling(R, N1L1620, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTNL/N1/L16to20/');
%N_fold_sampling(R, N1L1115, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTNL/N1/L11to15/');
%N_fold_sampling(R, N1L610, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTNL/N1/L6to10/');
%N_fold_sampling(R, N1L15, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTNL/N1/L1to5/');


N3T910=sparse(N3upTc9to10_rc(:,1),N3upTc9to10_rc(:,2),1,nchem,nprot);
N3T89=sparse(N3upTc8to9_rc(:,1),N3upTc8to9_rc(:,2),1,nchem,nprot);
N3T78=sparse(N3upTc7to8_rc(:,1),N3upTc7to8_rc(:,2),1,nchem,nprot);
N3T67=sparse(N3upTc6to7_rc(:,1),N3upTc6to7_rc(:,2),1,nchem,nprot);
N3T56=sparse(N3upTc5to6_rc(:,1),N3upTc5to6_rc(:,2),1,nchem,nprot);
N_fold_sampling(R, N3T910, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTMaxTc/N3/Tc9_10/');
N_fold_sampling(R, N3T89, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTMaxTc/N3/Tc8_9/');
N_fold_sampling(R, N3T78, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTMaxTc/N3/Tc7_8/');
N_fold_sampling(R, N3T67, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTMaxTc/N3/Tc6_7/');
N_fold_sampling(R, N3T56, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTMaxTc/N3/Tc5_6/');


N2T910=sparse(N2Tc9to10_rc(:,1),N2Tc9to10_rc(:,2),1,nchem,nprot);
N2T89=sparse(N2Tc8to9_rc(:,1),N2Tc8to9_rc(:,2),1,nchem,nprot);
N2T78=sparse(N2Tc7to8_rc(:,1),N2Tc7to8_rc(:,2),1,nchem,nprot);
N2T67=sparse(N2Tc6to7_rc(:,1),N2Tc6to7_rc(:,2),1,nchem,nprot);
N2T56=sparse(N2Tc5to6_rc(:,1),N2Tc5to6_rc(:,2),1,nchem,nprot);
N_fold_sampling(R, N2T910, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTMaxTc/N2/Tc9_10/');
N_fold_sampling(R, N2T89, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTMaxTc/N2/Tc8_9/');
N_fold_sampling(R, N2T78, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTMaxTc/N2/Tc7_8/');
N_fold_sampling(R, N2T67, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTMaxTc/N2/Tc6_7/');
N_fold_sampling(R, N2T56, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTMaxTc/N2/Tc5_6/');


N1T910=sparse(N1Tc9to10_rc(:,1),N1Tc9to10_rc(:,2),1,nchem,nprot);
N1T89=sparse(N1Tc8to9_rc(:,1),N1Tc8to9_rc(:,2),1,nchem,nprot);
N1T78=sparse(N1Tc7to8_rc(:,1),N1Tc7to8_rc(:,2),1,nchem,nprot);
N1T67=sparse(N1Tc6to7_rc(:,1),N1Tc6to7_rc(:,2),1,nchem,nprot);
N1T56=sparse(N1Tc5to6_rc(:,1),N1Tc5to6_rc(:,2),1,nchem,nprot);
N_fold_sampling(R, N1T910, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTMaxTc/N1/Tc9_10/');
N_fold_sampling(R, N1T89, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTMaxTc/N1/Tc8_9/');
N_fold_sampling(R, N1T78, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTMaxTc/N1/Tc7_8/');
N_fold_sampling(R, N1T67, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTMaxTc/N1/Tc6_7/');
N_fold_sampling(R, N1T56, 10, '/scratch/hansaim.lim/wizan/wiZAN/ZINC_data/chem_prot/NTMaxTc/N1/Tc5_6/');

end
