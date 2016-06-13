#!/usr/bin/python
import sys

ambig_file='./ZC_ambiguous_pairs.tsv'
ambig_ikey_uni=[]	#list of tuples (ikey, uniprot)
for line in open(ambig_file, "r").xreadlines():
	z=line.strip().split("\t")
	ikey=str(z[0])
	uni=str(z[1])
	pair=(ikey,uni)
	ambig_ikey_uni.append(pair)

infile='./ZC_ambiguous_pairs_assayInfo.tsv'
ikey_uniprot_ic50={}
for line in open(infile,"r").xreadlines():
        z=line.strip().split("\t")
	ikey=str(z[0])
	uni=str(z[1])
	ic50=float(z[2])
	ikey_uni=(ikey,uni)
	try:
		ikey_uniprot_ic50[ikey_uni].append(ic50)
	except:
		ikey_uniprot_ic50[ikey_uni]=[]
		ikey_uniprot_ic50[ikey_uni].append(ic50)

print "Standard_InChI_Key, UniProt_Accession, Ambiguity, Count_below10uM, Count_over10uM, max_IC50, min_IC50"
for key in ambig_ikey_uni:
	values=ikey_uniprot_ic50[key]
	ikey=str(key[0])
	uniprot=str(key[1])
	count_in=0	#number of IC50 values 10 micromolar (10000 nanomolar) or lower
	count_out=0
	ambiguity=None
	ic50_max=max(values)
	ic50_min=min(values)
	for val in values:
		if val > 10000: #values in nanomolar
			count_out += 1 #not good binding
		else:
			count_in += 1 #binding
	if count_in == 0:
		ambiguity='Inactive_all'
	elif count_out == 0:
		ambiguity='Active_all'
	else:
		ic50_out_near_threshold=sorted(values)[-count_out]	#nearest to 10 micromolar, but over the threshold
		ic50_in_near_threshold=sorted(values)[count_in - 1]		#nearest to 10 micromolar, but below the threshold
		if count_in > count_out:
			if ic50_max < 20000:
				#the highest IC50 is below 20 uM
				ambiguity='Active_tentative'
			else:
				ambiguity='Ambiguous_MaxIC50_too_high'
		else:
			ambiguity='Ambiguous_byMajorityVote'
	print "%s, %s, %s, %d, %d, %s, %s" % (ikey,uniprot,ambiguity,count_in,count_out,str(ic50_max),str(ic50_min))
