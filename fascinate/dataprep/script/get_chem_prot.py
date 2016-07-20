#!/usr/bin/python
import MySQLdb as db

chemprot_file='./combined_chem_prot_source.tsv'

con = db.connect('localhost', 'hlim', 'rhkdlf2043', 'hetio');
cur=con.cursor()
chem_qry="SELECT chemical_index, InChIKey, canonical_SMILES FROM chemical"
cur.execute(chem_qry)
chemline=cur.fetchall()
chem_by_ikey={}
for chem in chemline:
	ind=str(chem[0])
	ikey=str(chem[1])
	smi=str(chem[2])
	try:
		chem_by_ikey[ikey]
	except:
		chem_by_ikey[ikey]=[ind, smi]
prot_qry="SELECT gene_index, HGNC, gene_symbol, uniprot_accession, sequence_source, sequence FROM gene"
cur.execute(prot_qry)
protline=cur.fetchall()
prot_by_uniprot={}
for prot in protline:
	ind=str(prot[0])
	hgnc=str(prot[1])
	symbol=str(prot[2])
	uniprot=str(prot[3])
	src=str(prot[4])
	seq=str(prot[5])
	try:
		prot_by_uniprot[uniprot]
	except:
		prot_by_uniprot[uniprot]=[ind,hgnc,symbol,src,seq]

nf=open("./chem_prot_Not_Found.txt","w")
cpind=open("./chem_prot_byIndex.tsv","w")
for line in open(chemprot_file,"r").xreadlines():
	line=line.strip().split("\t")
	ikey=str(line[0])
	uniprot=str(line[1])
	cpsource=str(line[2])
	try:
		chem=chem_by_ikey[ikey]
		prot=prot_by_uniprot[uniprot]
		chem_index=chem[0]
		prot_index=prot[0]
		cpind.write(chem_index+"\t"+prot_index+"\n")
	except:
		nf.write(ikey+"\t"+uniprot+"\t"+cpsource+" Not Found by InChIKey and UniProt Accession\n")
		continue
con.close()	
