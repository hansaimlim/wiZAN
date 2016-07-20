#!/usr/bin/python

from Bio import SeqIO
import MySQLdb as db
import sys

con = db.connect('localhost', 'root', 'rhkdlf2043', 'hetio');
cur=con.cursor()
qry="SELECT uniprot_accession FROM gene"
cur.execute(qry)
result = cur.fetchall()
accs=[]
for s in result:
	accs.append(str(s[0]))


#sprots=[]	#list of gene symbols in uniprot_sprot
accs_nosprot=[]
count_sprot=0
handle=open("./uniprot_sprot.fasta","rU")
for record in SeqIO.parse(handle, "fasta"):
	ids=record.id.strip().split("|")
	acc=str(ids[1])
#	sprots.append(symbol)
	if acc in accs:
		#sequence is reviewed: source=Swiss-Prot
		qry1="UPDATE gene SET sequence = '%s', sequence_source='Swiss-Prot' WHERE uniprot_accession = '%s'" % (str(record.seq), acc)
		cur.execute(qry1)
		con.commit()
		count_sprot=count_sprot+1
	else:
		accs_nosprot.append(acc)
handle.close()

count_trembl=0
count_notfound=0
accs_notfound=[]
handle2=open("./uniprot_trembl.fasta","rU")
for record2 in SeqIO.parse(handle2, "fasta"):
	ids=record2.id.strip().split("|")
	acc=str(ids[1])
	if acc in accs_nosprot:
		qry2="UPDATE gene SET sequence = '%s', sequence_source='TrEMBL' WHERE uniprot_accession = '%s'" % (str(record2.seq), acc)
		cur.execute(qry2)
		con.commit()
		count_trembl=count_trembl+1
	else:
		accs_notfound.append(acc)
		count_notfound=count_notfound+1
total_found=count_sprot+count_trembl
print "sprot=%d, trembl=%d, notfound=%d, total_found=%d" % (count_sprot, count_trembl, count_notfound, total_found)
print accs_notfound

#for line in open('./diseases_sorted.txt',"r").xreadlines():
#	z=line.strip().split("\t")
#	doid=int(z[0])
#	name=str(db.escape_string(z[1])) or None
#	patho=str(db.escape_string(z[2])) or None
#	
#	qry="INSERT INTO disease (DOID, disease_name, pathophysiology) \
#	VALUES ('%d', '%s', '%s')" % (doid, name, patho)
#	try:
#		cur.execute(qry) 
#		con.commit()
#	except:
#		print z
#		con.rollback()
con.close()
