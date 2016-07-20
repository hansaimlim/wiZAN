#!/usr/bin/python
import MySQLdb as db


con = db.connect('localhost', 'hlim', 'rhkdlf2043', 'hetio');
cur=con.cursor()

ikeyfile='./ikey_notfound.txt'
unifile='./uniprot_notfound.txt'
for line in open(ikeyfile, "r").xreadlines():
	line=line.strip().split("\t")
	ikey=str(line[0])
	qry="SELECT cstr.standard_inchi_key, cstr.canonical_smiles, mdic.pref_name, mdic.chembl_id, crec.compound_name\
	     FROM chembl_20.compound_structures cstr INNER JOIN chembl_20.molecule_dictionary mdic ON cstr.molregno=mdic.molregno \
	     INNER JOIN chembl_20.compound_records crec ON cstr.molregno=crec.molregno\
	     WHERE cstr.standard_inchi_key = '%s'"%(ikey)
	try:
		cur.execute(qry)
		chemdata=cur.fetchone()
		ikey=str(chemdata[0])
		smi=str(chemdata[1])
		prefname=str(chemdata[2])
		chembl=str(chemdata[3])
		compname=str(chemdata[4])
		
		print "%s\t%s\t%s\t%s\t%s"%(ikey,smi,prefname,chembl,compname)
	except:
		None
con.close()	
