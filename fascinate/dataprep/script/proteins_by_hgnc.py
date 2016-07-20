#!/usr/bin/python

import httplib2 as http
import json
from Bio import SeqIO
import MySQLdb as db
import sys

con = db.connect('localhost', 'root', 'rhkdlf2043', 'hetio');
cur=con.cursor()
qry="SELECT HGNC FROM gene WHERE sequence IS NULL"
cur.execute(qry)
result = cur.fetchall()
hgnc_noseq=[]
for s in result:
        hgnc_noseq.append(int(s[0]))

try:
 from urlparse import urlparse
except ImportError:
 from urllib.parse import urlparse

accs=[]
for hgnc in hgnc_noseq:
	headers = {
	 'Accept': 'application/json',
	}

	uri = 'http://rest.genenames.org'
	method = 'GET'
	path = '/fetch/hgnc_id/'+str(hgnc)
	target = urlparse(uri+path)
	body = ''
	h = http.Http()
	response, content = h.request(
	 target.geturl(),
	 method,
	 body,
	 headers)

	if response['status'] == '200':
	 # assume that content is a json reply
	 # parse content with the json module 
		data = json.loads(content)
		try:
			acc=str(data['response']['docs'][0]['uniprot_ids'])
			accs.append(acc)

		except:
			print "Uniprot ID not found for HGNC: %d" % (hgnc)
	else:
		print 'Error detected: ' + response['status']

handle=open("./uniprot_sprot.fasta","rU")
for record in SeqIO.parse(handle, "fasta"):
	ids=record.id.strip().split("|")
	acc=str(ids[1])
#       sprots.append(symbol)
	if acc in accs:
		#sequence is reviewed: source=Swiss-Prot
		print "updating sequence info for %s" % (acc)
		qry1="UPDATE gene SET sequence = '%s', sequence_source='Swiss-Prot' WHERE uniprot_accession = '%s'" % (str(record.seq), acc)
		try:
			cur.execute(qry1)
			con.commit()
		except:
			print "sequence not found for: %s" % (acc)
handle.close()
