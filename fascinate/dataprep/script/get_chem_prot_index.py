#!/usr/bin/python
import MySQLdb as db

con = db.connect('localhost', 'hlim', 'w31c0m3', 'hetio');
cur=con.cursor()

def get_chemical_index_by_InChIKey(ikey):
	#ind is None if chemical not found
	qry="SELECT chemical_index FROM chemical WHERE InChIKey=%s"
	cur.execute(qry,(ikey))
	ind=cur.fetchone()
	return ind
def get_chemical_index_by_altInChIKey(ikey):
	#ind is None if chemical not found
	qry="SELECT chemical_index FROM chemical WHERE Alternate_InChIKey=%s"
	cur.execute(qry,(ikey))
	ind=cur.fetchone()
	return ind
def get_gene_index_by_UniProt(uni):
	qry="SELECT gene_index FROM gene WHERE uniprot_accession=%s"
	cur.execute(qry,(uni))
	ind=cur.fetchone()
	return ind
def insert_chem_prot(chemind,protind,source):
	chemind=str(chemind)
	protind=str(protind)
	source=str(source)
	qry="INSERT INTO chemical_gene(chemical_index,gene_index,source) VALUES (%s, %s, %s)"
	try:
		cur.execute(qry,(chemind,protind,source))
		con.commit()
	except:
		con.rollback()
def check_chem_prot(chemind,protind):
	chemind=str(chemind)
	protind=str(protind)
	qry="SELECT chemical_index, gene_index FROM chemical_gene WHERE chemical_index = %s AND gene_index = %s"
	cur.execute(qry,(chemind,protind))
	dat=cur.fetchall()
	return dat

S=open('../combined_chem_prot_index_add.csv',"w")
NF=open('../chemical_ikey_notfound_from_chemprot.txt',"w")
for line in open('../combined_chem_prot_source.tsv',"r").xreadlines():
	line=line.strip().split("\t")
	ikey=str(line[0])
	uni=str(line[1])
	source=str(line[2])
	chem=get_chemical_index_by_InChIKey(ikey)
	prot=get_gene_index_by_UniProt(uni)
	if check_chem_prot(chem,prot) is not None:
		#the pair already exists
		continue

	if chem is None:
		chem=get_chemical_index_by_altInChIKey(ikey)

	if chem is not None:
		if prot is not None:
			chemind=str(chem[0])
			protind=str(prot[0])
			#both chemical and protein exist in DB
			S.write(chemind+", "+protind+"\n")
			insert_chem_prot(chemind,protind,source)
		else:
			None
#			print "%s\tProtein Not found"%(uni)
	else:
		NF.write(ikey+"\n")
con.close()	
