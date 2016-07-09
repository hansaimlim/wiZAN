#!/usr/bin/python
import itertools

clu_file='./mcl_out_hierachy.tsv'
ccsim_file='./FDA_DrugNames_tanimoto_cosine.tsv'
skipped_edges=0
cluster_drug_pairs=[]
for line in open(clu_file,"r").xreadlines():
	r=line.strip().split("\t")
	for pair in itertools.combinations(r,2):
		cluster_drug_pairs.append(pair)

for line in open(ccsim_file,"r").xreadlines():
	r=line.strip().split("\t")
	drug1=str(r[0])
	drug2=str(r[1])
	tani =float(r[2])
	coss =float(r[3])
	pair1=(drug1,drug2)
	pair2=(drug2,drug1)
	
	if pair1 in cluster_drug_pairs:
		print "%s\t%s\t%s\t%s" %(drug1,drug2,tani,coss)
	elif pair2 in cluster_drug_pairs:
		print "%s\t%s\t%s\t%s" %(drug1,drug2,tani,coss)
	else:
		skipped_edges += 1
print "%d edges not found in cluster" % (skipped_edges)
