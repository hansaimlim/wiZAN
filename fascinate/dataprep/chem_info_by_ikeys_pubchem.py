#!/usr/bin/python

import httplib2 as http
import json
#import MySQLdb as db
import sys
try:
 from urlparse import urlparse
except ImportError:
 from urllib.parse import urlparse

#con = db.connect('localhost', 'root', 'rhkdlf2043', 'hetio');
#cur=con.cursor()
#qry="SELECT HGNC FROM gene WHERE sequence IS NULL"
#cur.execute(qry)
#result = cur.fetchall()
#hgnc_noseq=[]
#for s in result:
#        hgnc_noseq.append(int(s[0]))

###
### List of functions ###
# get_registrynumber_by_InChIKey(inchikey)
# get_canonicalsmiles_by_InChIKey(inchikey)
# get_CID_by_InChIKey(inchikey)
# get_synonym_by_InChIKey(inchikey)

# get_CID_by_CTD(meshID)
# get CID_by_name(chemical_name) 
# get_CID_by_CAS(CAS_number)
#
# get_InChIKey_by_CID(CID)
# get_canonicalsmiles_by_CID(CID)

def get_registrynumber_by_InChIKey(ikey):
	rn=None
	headers = {
	 'Accept': 'application/json',
	}
	method = 'GET'
	pre_uri='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/'
	suf_uri='/xrefs/rn/json'
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
			rn=str(data['InformationList']['Information'][0]['RN'][0])
		except:
			print "CAS RN not found from PubChem for: %s" % (ikey)
	return rn
def get_canonicalsmiles_by_InChIKey(ikey):
	smi=None
	headers = {
	 'Accept': 'application/json',
	}
	method = 'GET'
	pre_uri='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/'
	suf_uri='/property/canonicalsmiles/json'
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
			smi=str(data['PropertyTable']['Properties'][0]['CanonicalSMILES'])
		except:
			print "CanonicalSMILES not found from PubChem for: %s" % (ikey)
	return smi
def get_synonym_by_InChIKey(ikey):
	syn=None
	headers = {
	 'Accept': 'application/json',
	}
	method = 'GET'
	pre_uri='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/'
	suf_uri='/synonyms/json'
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
			syn=str(data['InformationList']['Information'][0]['Synonym'][0])
		except:
			print "Synonym not found from PubChem for: %s" % (ikey)
	return syn
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
			print "CID not found from PubChem for: %s" % (ikey)
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
def get_CID_by_name(name):
	cid=None
	headers = {
	 'Accept': 'application/json',
	}
	method = 'GET'
	pre_uri='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
	suf_uri='/cids/json'
	path = pre_uri + name + suf_uri
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
			print "CID not found from PubChem for: %s" % (name)
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


# for each inchikey, needed [synonym, CID, CAS RN, CanonicalSMILES]
chem_info=[]
infile=sys.argv[1]
outfile=sys.argv[2]
for line in open(infile,"r").xreadlines():
	line = line.strip().split("\t")
	ikey = str(line[0])
	cid=None
	syn=None
	rn=None
	smi=None
	try:
		cid=get_CID_by_InChIKey(ikey)
	except:
		cid=None
	try:
		smi=get_canonicalsmiles_by_InChIKey(ikey)
	except:
		smi=None
	if cid is not None:
	#cid is required
		if smi is not None:
		#canonical smiles is required
			try:
				syn=get_synonym_by_InChIKey(ikey)
			except:
				syn=None
			try:
				rn=get_registrynumber_by_InChIKey(ikey)
			except:
				rn=None
			cheminfo=[ikey, syn, cid, rn, smi]
			chem_info.append(cheminfo)
		

s=open(outfile,"w")
for cheminfo in chem_info:
	ikey=str(cheminfo[0])
	syn=str(cheminfo[1])
	cid=str(cheminfo[2])
	rn=str(cheminfo[3])
	smi=str(cheminfo[4])
	print "%s\t%s\t%s\t%s\t%s" % (ikey, syn, cid, rn, smi)
			
#handle=open("./uniprot_sprot.fasta","rU")
#for record in SeqIO.parse(handle, "fasta"):
#	ids=record.id.strip().split("|")
#	acc=str(ids[1])
#       sprots.append(symbol)
#	if acc in accs:
		#sequence is reviewed: source=Swiss-Prot
#		print "updating sequence info for %s" % (acc)
#		qry1="UPDATE gene SET sequence = '%s', sequence_source='Swiss-Prot' WHERE uniprot_accession = '%s'" % (str(record.seq), acc)
#		try:
#			cur.execute(qry1)
#			con.commit()
#		except:
#			print "sequence not found for: %s" % (acc)
#handle.close()
