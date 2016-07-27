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
cdfile='../output/more_chemicals.tsv'

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
def get_CAS_by_CID(cid):
        cas=None
        headers = {
         'Accept': 'application/json',
        }
        method = 'GET'
        pre_uri='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'
        suf_uri='/xrefs/rn/json'
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
def update_chem_CAS(chemind,cas):
        chemind=str(chemind)
        cas=str(cas)
        qry="UPDATE chemical SET CAS = %s WHERE chemical_index = %s"
        try:
                cur.execute(qry,(cid,chemind,))
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
def get_chemical_noCAS():
	chems=[]
	qry="SELECT chemical_index, InChIKey, PubChem_CID FROM chemical WHERE CAS IS NULL"
	cur.execute(qry)
	dat=cur.fetchall()
	for row in dat:
		ind=str(row[0])
		ikey=str(row[1])
		if row[2] is None:
			cid=None
		else:
			cid=str(row[2])
		r=(ind,ikey,cid)
		chems.append(r)
	return chems
def insert_chemical(ikey,altikey,chemname,cid,cas,chembl,altid,smiles):
        qry="INSERT INTO chemical(InChIKey,Alternate_InChIKey,chemical_name,PubChem_CID,CAS,ChEMBL_id,Alternate_id,canonical_SMILES)\
		 VALUES (%s,%s,%s,%s,%s,%s,%s,%s)"
        if get_chemical_index_by_InChIKey(ikey) is None:
                try:
                        cur.execute(qry,(ikey,altikey,chemname,cid,cas,chembl,altid,smiles,))
                        con.commit()
                except:
                        con.rollback()
def get_smiles_from_chembl(ikey):
        qry="SELECT cstr.canonical_smiles\
        FROM chembl_20.compound_structures cstr\
        INNER JOIN chembl_20.molecule_dictionary mdic ON cstr.molregno=mdic.molregno\
        WHERE cstr.standard_inchi_key=%s"
        cur.execute(qry,(ikey,))
        dat=cur.fetchone()
        if dat is None:
                return None
        else:
                return str(dat[0])
def get_synonym_from_chembl(ikey):
        qry="SELECT mdic.pref_name\
        FROM chembl_20.compound_structures cstr\
        INNER JOIN chembl_20.molecule_dictionary mdic ON cstr.molregno=mdic.molregno\
        WHERE cstr.standard_inchi_key=%s"
        cur.execute(qry,(ikey,))
        dat=cur.fetchone()
        if dat is None:
                return None
        else:
                return str(dat[0])
def get_chembl_id_from_chembl(ikey):
        qry="SELECT mdic.chembl_id\
        FROM chembl_20.compound_structures cstr\
        INNER JOIN chembl_20.molecule_dictionary mdic ON cstr.molregno=mdic.molregno\
        WHERE cstr.standard_inchi_key=%s"
        cur.execute(qry,(ikey,))
        dat=cur.fetchone()
        if dat is None:
                return None
        else:
                return str(dat[0])

def get_compound_name_from_chembl(ikey):
	qry="SELECT crec.compound_name FROM chembl_20.compound_structures cstr\
	INNER JOIN chembl_20.compound_records crec ON cstr.molregno=crec.molregno WHERE cstr.standard_inchi_key=%s"
        cur.execute(qry,(ikey,))
        dat=cur.fetchone()
        if dat is None:
                return None
        else:
                return str(dat[0])

linenum=1
for line in open(cdfile,"r").xreadlines():
	if linenum == 1:
		linenum+=1
		continue
	else:
		linenum+=1
	line=line.strip().split("\t")
	ikey=str(line[0])
	altikey=str(line[1])
	syn=str(line[2])
	cid=str(line[3])
	cas=str(line[4])
	chembl=str(line[5])
	altid=str(line[6])
	smi=str(line[7])
	
	if (altikey=='') or (altikey is None):
		altikey=None
	if (syn=='') or (syn is None):
		syn=get_synonym_from_chembl(ikey)
	if (cid=='') or (cid is None):
		cid=None
	if (cas=='') or (cas is None):
		if cid is not None:
			cas=get_CAS_by_CID(cid)
		else:
			cas=None
	if (chembl=='') or (chembl is None):
		chembl=get_chembl_id_from_chembl(ikey)
	if (altid=='') or (altid is None):
		altid=None	
	if syn==chembl:
		syn=get_compound_name_from_chembl(ikey)
	insert_chemical(ikey,altikey,chemname,cid,cas,chembl,altid,smiles)
	if (linenum % 2000) == 0:
		print "%d lines processed"%(linenum)

print "Processing CAS update"
chems=get_chemical_noCAS()
cid_count=0
cas_count=0
for chem in chems:
        ind=chem[0]
        ikey=chem[1]
        cid=chem[2]
	if cid is None:
		cid=get_CID_by_InChIKey(ikey)
		if cid is not None:
			#new CID found, update
			update_chem_CID(ind,cid)
			cid_count+=1
	if cid is None:
		None
	else:
		cas=get_CAS_by_CID(cid)
		if cas is not None:
			#CAS found. update
			update_chem_CAS(ind,cas)
			cas_count+=1
con.close()	
print "%d CID updated"%(cid_count)
print "%d CAS updated"%(cas_count)
