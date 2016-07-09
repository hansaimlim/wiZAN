#!/usr/bin/python
import itertools
#from collections import Counter

first_clu_file='./chem_chem/mcl_out_hierachy_first.tsv'
gpcr_file='./output/FDAdrug_GPCR_NormScore_Knownassociation.tsv'
kinase_file='./output/FDAdrug_Kinase_NormScore_Knownassociation.tsv'
gpcr_out='./output/first_cluster_GPCRs.tsv'
kinase_out='./output/first_cluster_Kinases.tsv'
not_found='./output/first_cluster_drugs_noGPCR_noKinase.txt'
drug_gpcr={}
drug_kinase={}
gpcr_count=0
kinase_count=0


for line in open(gpcr_file,"r").xreadlines():
	r=line.strip().split("\t")
	Drugname=str(r[0])
	drugname=Drugname.lower()
	target=str(r[1])
	score=float(r[2])
	known=int(r[3])
	if known == 1:
		try:
			drug_gpcr[drugname].append(target)
		except:
			drug_gpcr[drugname]=[target]

for line in open(kinase_file,"r").xreadlines():
	r=line.strip().split("\t")
	Drugname=str(r[0])
	drugname=Drugname.lower()
	target=str(r[1])
	score=float(r[2])
	known=int(r[3])
	if known == 1:
		try:
			drug_kinase[drugname].append(target)
		except:
			drug_kinase[drugname]=[target]

GPCR=open(gpcr_out,"w")
Kinase=open(kinase_out,"w")
NF=open(not_found,"w")
sep=";"
for line in open(first_clu_file,"r").xreadlines():
	drugs=line.strip().split("\t")
	for dru in drugs:
		drug=dru.lower()
		try:
			gpcrs=drug_gpcr[drug]
			gpcr_count=gpcr_count+len(gpcrs)
			gpcr_targets=""
			for gp in gpcrs:
				gpcr_targets=gpcr_targets+", "+gp
			GPCR.write(str(drug)+"\t"+str(gpcr_targets)+"\n")
		except:
			NF.write(str(drug)+" does not have known GPCR targets\n")
		try:
			kinases=drug_kinase[drug]
			kinase_count=kinase_count+len(kinases)
			kinase_targets=""
			for ki in kinases:
				kinase_targets=kinase_targets+", "+ki
			Kinase.write(str(drug)+"\t"+str(kinase_targets)+"\n")
		except:
			NF.write(str(drug)+" does not have known Kinase targets\n")
				
print "%d GPCRs found"%(gpcr_count)
print "%d Kinases found"%(kinase_count)
