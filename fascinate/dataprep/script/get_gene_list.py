#!/usr/bin/python

import MySQLdb as db
import sys

con = db.connect('localhost', 'root', 'rhkdlf2043', 'hetio');
cur=con.cursor()
qry="SELECT gene_index, ensembl, uniprot_accession, HGNC FROM gene WHERE sequence IS NOT NULL"
cur.execute(qry)
result = cur.fetchall()
for s in result:
	print "%s\t%s\t%s\t%s" % (str(s[0]),str(s[1]),str(s[2]),str(s[3]))
