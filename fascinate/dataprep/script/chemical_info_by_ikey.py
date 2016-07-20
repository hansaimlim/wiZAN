#!/usr/bin/python

import httplib2 as http
import json
import MySQLdb as db
import sys

#con = db.connect('localhost', 'root', 'rhkdlf2043', 'hetio');
#cur=con.cursor()
#qry="SELECT HGNC FROM gene WHERE sequence IS NULL"
#cur.execute(qry)
#result = cur.fetchall()
#hgnc_noseq=[]
#for s in result:
#        hgnc_noseq.append(int(s[0]))
def get_CID_by_InChIKey(ikey):
	cid=None
	headers = {
	 'Accept': 'application/json',
	}
	method = 'GET'
	pre_uri='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/'
	suf_uri='/cids/json'
	path = pre_uri + ikey + suf_uri
	target = urlparse(path)
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
			cid=str(data['IdentifierList']['CID'][0])
		except:
			print "CID not found from PubChem for: %s" % (cas)
	return cid

def get_CID_by_CTD(meshID):
	cid=None
	headers = {
	 'Accept': 'application/json',
	}
	method = 'GET'
	pre_uri='https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sourceid/Comparative%20Toxicogenomics%20Database/'
	suf_uri='/cids/json'
	path = pre_uri + meshID + suf_uri
	target = urlparse(path)
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
			cid=str(data['InformationList']['Information'][0]['CID'][0])
		except:
			print "CID not found from PubChem for: %s" % (meshID)
	return cid

def get_CID_by_CAS(cas):
	cid=None
	headers = {
	 'Accept': 'application/json',
	}
	method = 'GET'
	pre_uri='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/registryid/'
	suf_uri='/cids/json'
	path = pre_uri + cas + suf_uri
	target = urlparse(path)
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
			cid=str(data['IdentifierList']['CID'][0])
		except:
			print "CID not found from PubChem for: %s" % (cas)
	return cid

def get_InChIKey_by_CID(cid):
	ikey=None
	headers = {
	 'Accept': 'application/json',
	}
	method = 'GET'
	pre_uri='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'
	suf_uri='/property/InChIKey/json'
	path = pre_uri + str(cid) + suf_uri
	target = urlparse(path)
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
			ikey=str(data['PropertyTable']['Properties'][0]['InChIKey'])
		except:
			print "InChIKey not found from PubChem for: %s" % (str(cid))
	return ikey
def get_canonicalsmiles_by_CID(cid):
	smiles=None
	headers = {
	 'Accept': 'application/json',
	}
	method = 'GET'
	pre_uri='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'
	suf_uri='/property/CanonicalSMILES/json'
	path = pre_uri + str(cid) + suf_uri
	target = urlparse(path)
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
			smiles=str(data['PropertyTable']['Properties'][0]['CanonicalSMILES'])
		except:
			print "Canonical SMILES not found from PubChem for: %s" % (str(cid))
	return smiles


try:
 from urlparse import urlparse
except ImportError:
 from urllib.parse import urlparse

#need to know [cid, inchikey, chemical name, cas number, canonical smiles]
for line in open('./uniq_ikeys.txt',"r").xreadlines():
	line=line.strip().split("\t")
	ikey=str(line[0])
	try:
		cid=get_CID_by_InChIKey(ikey)
	except:
		cid=None
	print ikey

#s=open('./CTD_drugs_smiles.tsv',"w")
#synout=open('./CTD_drugs_synonyms.tsv',"w")
#for ctd in ctd_ikey_smile:
#	val=ctd_ikey_smile[ctd]
#	ikey=val[0]
#	smiles=val[1]
#	if ikey_count[ikey]==1:
#		s.write(ctd+"\t"+ikey+"\t"+smiles+"\n")
#	else:
#		syn=''
#		syns=ikey_synonyms[ikey]
#		for i in range(len(syns)):
#			syn=syn+syns[i]+"\t"
#		synout.write(syn+"\n")
			
