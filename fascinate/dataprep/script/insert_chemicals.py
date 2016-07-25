#!/usr/bin/python
import MySQLdb as db
import httplib2 as http
import json
import sys
try:
 from urlparse import urlparse
except ImportError:
 from urllib.parse import urlparse

con = db.connect('localhost', 'hlim', 'w31c0m3', 'hetio');
cur=con.cursor()

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
			None
                       #print "Synonym not found from PubChem for: %s" % (ikey)
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
                        None
                       # print "CID not found from PubChem for: %s" % (ikey)
        return cid
def get_CAS_by_InChIKey(ikey):
        cas=None
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
                        cas=str(data['InformationList']['Information'][0]['RN'][0])
                except:
                        None
                       # print "CanonicalSMILES not found from PubChem for: %s" % (ikey)
        return cas
def insert_chemical(ikey,altikey,chemname,cid,cas,chembl,altid,smiles):
        qry="INSERT INTO chemical(InChIKey,Alternate_InChIKey,chemical_name,PubChem_CID,CAS,ChEMBL_id,Alternate_id,canonical_SMILES)\
                 VALUES (%s,%s,%s,%s,%s,%s,%s,%s)"
        if get_chemical_index_by_InChIKey(ikey) is None:
                try:
                        cur.execute(qry,(ikey,altikey,chemname,cid,cas,chembl,altid,smiles))
                        con.commit()
                except:
                        con.rollback()


pubchemfile='../output/chem_info_from_PubChem.tsv'
smifile='../output/ikey_infoFound_smilesOnly.txt'
drugfile='../output/ikey_infoFound_drugname.txt'

for line in open(pubchemfile,"r").xreadlines():
	line=line.strip().split("\t")
	ikey=str(line[0])
	if line[1] is not None:
		syn=str(line[1])
	else:
		syn=None
	if line[2] is not None:
		cid=str(line[2])
	else:
		cid=None
	if line[3] is not None:
		cas=str(line[3])
	else:
		cas=None
	smi=str(line[4])
	altikey=None
	chembl=None
	altid=None
	insert_chemical(ikey,altikey,syn,cid,cas,chembl,altid,smi)

for  line in open(drugfile,"r").xreadlines():
	line=line.strip().split("\t")
	ikey=str(line[0])
	syn=str(line[1])
	cid=None
	try:
		cid=get_CID_by_InChIKey(ikey)
	except:
		cid=None
	cas=None
	try:
		cas=get_CAS_by_InChIkey(ikey)
	except:
		cas=None
	smi=str(line[2])
	altikey=None
	chembl=None
	altid=None
	insert_chemical(ikey,altikey,syn,cid,cas,chembl,altid,smi)
	
			
for  line in open(smifile,"r").xreadlines():
	line=line.strip().split("\t")
	ikey=str(line[0])
	syn=None
	try:
		syn=get_synonym_by_InChIKey(ikey)
	except:
		None
	cid=None
	try:
		cid=get_CID_by_InChIKey(ikey)
	except:
		cid=None
	cas=None
	try:
		cas=get_CAS_by_InChIkey(ikey)
	except:
		cas=None
	smi=str(line[1])
	altikey=None
	chembl=None
	altid=None
	insert_chemical(ikey,altikey,syn,cid,cas,chembl,altid,smi)

#chemfile='../chems_from_chembl.tsv'
#for line in open(chemfile, "r").xreadlines():
#        line=line.strip().split("\t")
#        ikey=str(line[0])
#        smi=str(line[1])
#	if line[2] is not None:
#	        prefname=str(line[2]) #often None
#	else:
#		prefname=None
#        chembl=str(line[3])
#        compname=str(line[4]) #long compound name
#	altikey=None
#	altid=None
#	insert_chemical(ikey,altikey,prefname,cid,cas,chembl,altid,smi)
#	
con.close()	
