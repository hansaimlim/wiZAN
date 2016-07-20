#!/usr/bin/python
import MySQLdb as db
import httplib2 as http
import json
import sys
try:
 from urlparse import urlparse
except ImportError:
 from urllib.parse import urlparse

con = db.connect('localhost', 'hlim', 'rhkdlf2043', 'hetio');
cur=con.cursor()

def get_disease_index_by_name(name):
	qry="SELECT disease_index FROM disease WHERE disease_name=%s"
	cur.execute(qry,(name))
	ind=cur.fetchone()
	return ind
def get_disease():
	qry="SELECT disease_index, disease_name FROM disease"
	cur.execute(qry)
	dat=cur.fetchall()
	name_index={}
	for row in dat:
		ind=str(row[0])
		name=str(row[1])
		name=name.lower()
		try:
			name_index[name]
		except:
			name_index[name]=ind
	return name_index
		

diseasename_index=get_disease()


for line in open('./therapeuticChemical_disease_CTD.txt',"r").xreadlines():
	line=line.strip().split("\t")
	diseasename_lower=str(line[0])
	chemicalname=str(line[1])
	chemId=str(line[2])
	cas=str(line[3])
	DiseaseName=str(line[4])
	disease_MESH=str(line[5])
	disind=diseasename_index[diseasename_lower]
	print disind

con.close()	
