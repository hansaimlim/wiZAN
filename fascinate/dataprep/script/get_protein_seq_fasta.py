#!/usr/bin/python

from Bio import SeqIO
import MySQLdb as db
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
con = db.connect('localhost', 'hlim', 'w31c0m3', 'hetio');
cur=con.cursor()
qry="SELECT gene_index, HGNC, gene_symbol, uniprot_accession, protein_name, chromosome, sequence FROM gene"
cur.execute(qry)
rows = cur.fetchall()
prots=[]
rownum=0
for row in rows:
	rownum+=1
	if rownum >5:
		break

	gi=str(row[0])
	hgnc="HGNC:"+str(row[1])
	symbol=row[2]
	acc=row[3]
	protname=row[4]
	chromosome=row[5]
	sequence=str(row[6])
	identifiers=gi+"|"+hgnc+"|"+symbol+"|"+acc
	desc=protname+"|"+chromosome
	tup=(gi,hgnc,symbol,acc,protname,chromosome,sequence)
	record=SeqRecord(Seq(sequence,IUPAC.protein),id=identifiers,description=desc,name=protname)
	print record
#	prots.append(tup)
	

#sprots=[]	#list of gene symbols in uniprot_sprot
con.close()
