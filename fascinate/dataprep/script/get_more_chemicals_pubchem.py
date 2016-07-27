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

def get_chemical_index_by_InChIKey(ikey):
        #ind is None if chemical not found
        qry="SELECT chemical_index FROM chemical WHERE InChIKey=%s"
        cur.execute(qry,(ikey,))
        ind=cur.fetchone()
        return ind

def get_chemicals_ikey():
        ikeys=[]
        qry="SELECT chemical_index, InChIKey FROM chemical"
        cur.execute(qry,)
        dat=cur.fetchall()
	for chem in dat:
		idx=str(dat[0])
		ikey=str(dat[1])
		ikeys.append(ikey)
        return ikeys
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
                return str(dat)
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
                return str(dat)
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
                return str(dat)
ikeys=get_chemicals_ikey()
addfile='../output/ikey_infoNotFound.txt'
count_new=0
smi_notfound=0
print "InChIKey\tAlt_InChIKey\tSynonym\tPubChem_CID\tCAS\tChEMBL\tAlt_ID\tSMILES"
for line in open(addfile,"r").xreadlines():
	line=line.strip().split("\t")
	ikey=str(line[0])
	if ikey in ikeys:
#		print "%s found in DB" % (ikey)
		continue
	else:
		count_new+=1
	try:
		smi=get_canonicalsmiles_by_InChIKey(ikey)
	except:
		smi=None
	if smi is None:
		try:
			smi=get_smiles_from_chembl(ikey)
			smi=str(smi)	#does not work if smi is None
		except:
			smi_notfound+=1
			continue
	try:
		syn=get_synonym_by_InChIKey(ikey)
	except:
		syn=get_synonym_from_chembl(ikey)
	try:
		cid=get_CID_by_InChIKey(ikey)
	except:
		cid=None
	try:
		cas=get_CAS_by_InChIkey(ikey)
	except:
		cas=None
	altikey=None
	chembl=get_chembl_id_from_chembl(ikey)
	altid=None
	if ikey is None:
		ikey=''
	if altikey is None:
		altikey=''
	if syn is None:
		syn=''
	if cid is None:
		cid=''
	if cas is None:
		cas=''
	if chembl is None:
		chembl=''
	if altid is None:
		altid=''
	if smi is None:
		smi=''
#	insert_chemical(ikey,altikey,syn,cid,cas,chembl,altid,smi)
	print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(str(ikey), str(altikey), str(syn), str(cid), str(cas), str(chembl), str(altid), str(smi))
#print "%d new chemicals; $d chemicals without smiles"%(count_new,smi_notfound)
con.close()	
