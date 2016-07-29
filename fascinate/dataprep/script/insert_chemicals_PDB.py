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
                        try:
                                cas=str(data['InformationList']['Information'][1]['RN'][0])
                        except:
                                try:
                                        cas=str(data['InformationList']['Information'][2]['RN'][0])
                                except:
                                        cas=None

        return cas
def insert_chemical(ikey,altikey,chemname,cid,cas,chembl,altid,smiles):
	val=(ikey,altikey,chemname,cid,cas,chembl,altid,smiles,)
        qry="INSERT INTO chemical(InChIKey,Alternate_InChIKey,chemical_name,PubChem_CID,CAS,ChEMBL_id,Alternate_id,canonical_SMILES)\
                 VALUES (%s,%s,%s,%s,%s,%s,%s,%s)"
        if get_chemical_index_by_InChIKey(ikey) is None:
                try:
                        cur.execute(qry,val)
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


def get_chemical_index_by_InChIKey(ikey):
        #ind is None if chemical not found
        qry="SELECT chemical_index FROM chemical WHERE InChIKey=%s"
        cur.execute(qry,(ikey,))
        ind=cur.fetchone()
        return ind


chemfile='../chems_from_PDB.tsv'
linenum=1
for line in open(chemfile, "r").xreadlines():
	if linenum==1:
		linenum+=1
		continue
        line=line.strip().split("\t")
        ikey=str(line[0])
        smi=str(line[1])
	try:
		if line[2] is not None:
			prefname=str(line[2]) #often None
		else:
			if line[3] is not None:
				prefname=str(line[3])
			else:
				prefname=None
	except:
		prefname=get_synonym_from_chembl(ikey)
        chembl=get_chembl_id_from_chembl(ikey)
        compname=None
	altikey=None
	altid=None
	insert_chemical(ikey,altikey,prefname,cid,cas,chembl,altid,smi)
con.close()
print "Inserting chemicals from PDB is done"
