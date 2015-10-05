#!/usr/bin/python
import sys
import MySQLdb
import csv

chem_index_name={}
chem_index_ikey={}
mapfile='./DrugBank_FDAdrugs_index_inchikey_name.tsv'
with open(mapfile, "rb") as mp:
	for line in mp:
		row=line.strip().split("\t")
		index=int(row[0])
		ikey=str(row[1])
		name=str(row[2])
		chem_index_name[index]=name
		chem_index_ikey[index]=ikey
infile='./chem_chem/FDAdrug_TanimotoSim_CosineSim.tsv'
with open(infile, "rb") as sim:
	for line in sim:
		row=line.strip().split("\t")
		index1=int(row[0])
		index2=int(row[1])
		tanimoto=float(row[2])
		cossim=float(row[3])

		drug1=str(chem_index_name[index1])
		drug2=str(chem_index_name[index2])
		ikey1=str(chem_index_ikey[index1])
		ikey2=str(chem_index_ikey[index2])
		output="%s\t%s\t%f\t%f" % (drug1, drug2, tanimoto, cossim)
		print output
