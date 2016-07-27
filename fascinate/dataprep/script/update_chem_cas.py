#!/usr/bin/python

#update second CAS for chemicals having 1 CAS RN
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
			None
                       # print "CID not found from PubChem for: %s" % (ikey)
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
			None
                       # print "CID not found from PubChem for: %s" % (meshID)
        return cid
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
			None
                       # print "CanonicalSMILES not found from PubChem for: %s" % (ikey)
        return smi
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
			try:
        	                cas=str(data['InformationList']['Information'][1]['RN'][0])
			except:
				try:
					cas=str(data['InformationList']['Information'][2]['RN'][0])
				except:
					cas=None
					
        return cas
def get_InChIKey_by_name(name):
        ikey=None
        headers = {
         'Accept': 'application/json',
        }
        method = 'GET'
        pre_uri='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
        suf_uri='/property/inchikey/json/'
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
                        ikey=str(data['PropertyTable']['Properties'][0]['InChIKey'])
                except:
			None
        return ikey
def get_InChIKey_by_CID(cid):
        ikey=None
        headers = {
         'Accept': 'application/json',
        }
        method = 'GET'
        pre_uri='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'
        suf_uri='/property/inchikey/json/'
        path = pre_uri + cid + suf_uri
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
			None
        return ikey
def get_alt_InChIKey_by_CID(cid):
        ikey=None
        headers = {
         'Accept': 'application/json',
        }
        method = 'GET'
        pre_uri='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'
        suf_uri='/property/inchikey/json/'
        path = pre_uri + cid + suf_uri
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
                        ikey=str(data['PropertyTable']['Properties'][1]['InChIKey'])
                except:
			None
        return ikey
def get_InChIKey_by_CAS(cas):
	cas=str(cas)
        ikey=None
        headers = {
         'Accept': 'application/json',
        }
        method = 'GET'
        pre_uri='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/rn/'
        suf_uri='/property/inchikey/json/'
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
                        ikey=str(data['PropertyTable']['Properties'][0]['InChIKey'])
                except:
			None
        return ikey
	

def get_chemical_index_by_InChIKey(ikey):
	#ind is None if chemical not found
	qry="SELECT chemical_index FROM chemical WHERE InChIKey=%s"
	cur.execute(qry,(ikey,))
	ind=cur.fetchone()
	return ind
def get_chemical_index_by_CID(cid):
	#ind is None if chemical not found
	qry="SELECT chemical_index FROM chemical WHERE PubChem_CID=%s"
	cur.execute(qry,(cid,))
	ind=cur.fetchone()
	return ind
def get_chemical_index_by_CAS(cas):
	#ind is None if chemical not found
	qry="SELECT chemical_index FROM chemical WHERE CAS=%s"
	cur.execute(qry,(cas,))
	ind=cur.fetchone()
	return ind
def get_disease_index_by_name(name):
	qry="SELECT disease_index FROM disease WHERE disease_name=%s"
	cur.execute(qry,(name,))
	ind=cur.fetchone()
	return ind
def insert_chem_disease(chemind,disind):
	chemind=str(chemind)
	disind=str(disind)
	qry="INSERT INTO chemical_disease(chemical_index,disease_index) VALUES (%s, %s)"
	try:
		cur.execute(qry,(chemind,disind,))
		con.commit()
	except:
		con.rollback()
def update_chem_altid(chemind,altid):
	chemind=str(chemind)
	altid=str(altid)
	qry="UPDATE chemical SET Alternate_id = %s WHERE chemical_index = %s"
	try:
		cur.execute(qry,(altid,chemind,))
		con.commit()
	except:
		con.rollback()
def update_chem_altikey(chemind,altikey):
	chemind=str(chemind)
	altikey=str(altikey)
	qry="UPDATE chemical SET Alternate_InChIKey = %s WHERE chemical_index = %s"
	try:
		cur.execute(qry,(altikey,chemind,))
		con.commit()
	except:
		con.rollback()
def update_chem_CID(chemind,cid):
	chemind=str(chemind)
	cid=str(cid)
	qry="UPDATE chemical SET PubChem_CID = %s WHERE chemical_index = %s"
	try:
		cur.execute(qry,(cid,chemind,))
		con.commit()
	except:
		con.rollback()
def update_chem_CAS(chemind,cas):
	chemind=str(chemind)
	cas=str(cas)
	qry="UPDATE chemical SET CAS = %s WHERE chemical_index = %s"
	try:
		cur.execute(qry,(cid,chemind,))
		con.commit()
	except:
		con.rollback()
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
def insert_chemical(ikey,altikey,chemname,cid,cas,chembl,altid,smiles):
        qry="INSERT INTO chemical(InChIKey,Alternate_InChIKey,chemical_name,PubChem_CID,CAS,ChEMBL_id,Alternate_id,canonical_SMILES)\
		 VALUES (%s,%s,%s,%s,%s,%s,%s,%s)"
        if get_chemical_index_by_InChIKey(ikey) is None:
                try:
                        cur.execute(qry,(ikey,altikey,chemname,cid,cas,chembl,altid,smiles,))
                        con.commit()
                except:
                        con.rollback()
get_chems_qry="SELECT chemical_index, InChIKey FROM chemical WHERE CAS IS NULL" 
cur.execute(get_chems_qry)
chems=cur.fetchall()
for chem in chems:
	chemind=str(chem[0])
	ikey=str(chem[1])
	cas=get_CAS_by_InChIKey(ikey)
	if cas is not None:
		update_chem_CAS(chemind,cas)

print "CAS update is done"

#	if altikey is None:
#		print "Alternate ikey is null"
#	
#	if altid is None:
#		print "Alternate id is null"
con.close()	
