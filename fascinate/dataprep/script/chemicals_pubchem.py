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

cname_cid_ikey={}
ctd_ikey_smile={}
ikey_count={}
ikey_synonyms={}
for line in open('./CTD_drug.txt',"r").xreadlines():
	z = line.strip().split("\t")
	cname = str(z[0])
	if len(z) == 3:
		#z[1]=MeSH; z[2]=CAS; three searches available
		ctd=str(z[1])
		cas=str(z[2])
		cid=None
		try:
			cid=get_CID_by_CTD(ctd)
		except:
			try:
				cid=get_CID_by_CAS(cas)
			except:
				print "Could not find CID from PubChem for %s;%s" % (ctd,cas)
		if cid is not None:
			try:
				ikey=get_InChIKey_by_CID(cid)
				smiles=get_canonicalsmiles_by_CID(cid)
				ctd_ikey_smile[ctd]=(ikey,smiles)
			except:
				print "Either InChIKey or SMILES (or both) is not available for CID: %s" % (cid)
			try:
				ikey_count[ikey]
			except KeyError:
				ikey_count[ikey]=1
			else:
				ikey_count[ikey] += 1

			if ikey_count[ikey]==1:
				#ikey appearing first time
				ikey_synonyms[ikey]=[ctd,cas]
			else:
				ikey_synonyms[ikey].append(ctd)
				ikey_synonyms[ikey].append(cas)
			
				
	elif len(z) == 2:
		#no CAS
		ctd=str(z[1])
		cid=None
		try:
			cid=get_CID_by_CTD(ctd)
		except:
			print "Could not find CID from PubChem for %s" % (ctd)
		if cid is not None:
			try:
				ikey=get_InChIKey_by_CID(cid)
				smiles=get_canonicalsmiles_by_CID(cid)
				ctd_ikey_smile[ctd]=(ikey,smiles)
			except:
				print "Either InChIKey or SMILES (or both) is not available for CID: %s" % (cid)
			try:
				ikey_count[ikey]
			except KeyError:
				ikey_count[ikey]=1
			else:
				ikey_count[ikey] += 1

			if ikey_count[ikey] == 1:
				#ikey appearing first time
				ikey_synonyms[ikey]=[ctd]
			else:
				ikey_synonyms[ikey].append(ctd)
	else:
		#cname only
		print "no other identifiers available than the name"

s=open('./CTD_drugs_smiles.tsv',"w")
synout=open('./CTD_drugs_synonyms.tsv',"w")
for ctd in ctd_ikey_smile:
	val=ctd_ikey_smile[ctd]
	ikey=val[0]
	smiles=val[1]
	if ikey_count[ikey]==1:
		s.write(ctd+"\t"+ikey+"\t"+smiles+"\n")
	else:
		syn=''
		syns=ikey_synonyms[ikey]
		for i in range(len(syns)):
			syn=syn+syns[i]+"\t"
		synout.write(syn+"\n")
			
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
