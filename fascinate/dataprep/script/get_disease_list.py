#!/usr/bin/python

import MySQLdb as db
import sys

con = db.connect('localhost', 'root', 'rhkdlf2043', 'hetio');
cur=con.cursor()
qry="SELECT disease_index, disease_name FROM disease"
cur.execute(qry)
result = cur.fetchall()
for s in result:
	print "%s\t%s" % (str(s[0]),str(s[1]))
