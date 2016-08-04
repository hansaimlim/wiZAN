parfor i=1:3
	rank=200;
	output_dir='/home/hlim/wizan/wiZAN/ZINC_data/TPR/';
	path_pref='/home/hlim/wizan/wiZAN/ZINC_data/chem_prot/';
	ntnl_suff='numtarget_numligand/';
	ntmaxtc_suff='numtarget_maxTc/';
	nt=['N1','N2','N3up'];
	csvsuff='.csv'
	ntnl=['L1to5','L6to10','L11to15','L16to20','L21more'];
	maxtc=['Tc0.49to0.6','Tc0.6to0.7','Tc0.7to0.8','Tc0.8to0.9','Tc0.9to1'];
	NT=nt(i);
	for j=1:5
		NTNL_train=[path_pref ntnl_suff 'train_' NT ntnl(j) csvsuff];
		NTNL_test=[path_pref ntnl_suff 'test_' NT ntnl(j)] csvsuff;
		NTNL_out=[output_dir NT ntnl(j) '_TPR.txt'];

		NTmaxTc_train=[path_pref ntmaxtc_suff 'train_' NT maxtc(j) csvsuff];
		NTmaxTc_test=[path_pref ntmaxtc_suff 'test_' NT maxtc(j) csvsuff];
		NTmaxTc_out=[output_dir NT maxtc(j) '_TPR.txt'];

		wiZAN_dual_csv(NTNL_train, NTNL_test, rank, NTNL_out);
		wiZAN_dual_csv(NTmaxTc_train, NTmaxTc_test, rank, NTmaxTc_out);
	end		
end		
