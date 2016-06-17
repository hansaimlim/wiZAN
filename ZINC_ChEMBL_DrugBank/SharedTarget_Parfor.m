function ShareStat=SharedTarget_Parfor()
%function to get chemical-chemical structural similarity vs. shared target ratio statistics
%each bin contains the shared-target ratio for the given chem-chem pairs within the similarity range
%each data point is for each chem-chem pair, NOT the individual chemical
clear
tic;
load /scratch/hansaim.lim/wizan/wiZAN/ZINC_ChEMBL_DrugBank/chem_chem/chem_chem_ZCD;
load /scratch/hansaim.lim/wizan/wiZAN/ZINC_ChEMBL_DrugBank/chem_prot/chem_prot_ZCD;
chem_chem=chem_chem_ZCD;
chem_prot=chem_prot_ZCD;
total_chem=0;total_target=0;num_chem_shared=0;num_target_shared=0;
end_list=[70500, 99700, 122000, 141000, 157600, 172600, 186400, 199338];

parpool(numel(end_list));

parfor i=numel(end_list)
 if i==1
  start_index=1;
 else
  start_index=end_list(i-1)+1;
 end
 for r=start_index:end_list(i)
  for c=r+1:end_list(i)
   if chem_chem(r,c)<0.5 % the chemical-chemical similarity condition
    total_chem=total_chem+1;
    cp1=chem_prot_ZCD(r,:);
    cp2=chem_prot_ZCD(c,:);
    cps=cp1+cp2;
    tot_tar=sum(cps>0);
    share_tar=sum(cps>1);
    total_target=total_target+tot_tar;
    if share_tar > 0
     num_chem_shared=num_chem_shared+1;
     num_target_shared=num_target_shared+share_tar;
    end
   end
  end
 end
end

ShareStat=[total_chem, num_chem_shared, total_target, num_target_shared];
end
