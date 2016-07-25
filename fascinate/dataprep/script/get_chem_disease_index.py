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

def get_chemical_index_by_InChIKey(ikey):
	#ind is None if chemical not found
	qry="SELECT chemical_index FROM chemical WHERE InChIKey=%s"
	cur.execute(qry,(ikey))
	ind=cur.fetchone()
	return ind
def get_chemical_index_by_CID(cid):
	#ind is None if chemical not found
	qry="SELECT chemical_index FROM chemical WHERE PubChem_CID=%s"
	cur.execute(qry,(cid))
	ind=cur.fetchone()
	return ind
def get_chemical_index_by_CAS(cas):
	#ind is None if chemical not found
	qry="SELECT chemical_index FROM chemical WHERE CAS=%s"
	cur.execute(qry,(cas))
	ind=cur.fetchone()
	return ind
def get_disease_index_by_name(name):
	qry="SELECT disease_index FROM disease WHERE disease_name=%s"
	cur.execute(qry,(name))
	ind=cur.fetchone()
	return ind
def insert_chem_disease(chemind,disind):
	chemind=str(chemind)
	disind=str(disind)
	qry_check="SELECT * FROM chemical_disease WHERE chemical_index = %s AND disease_index = %s"
	cur.execute(qry_check,(chemind,disind))
	check=cur.fetchone()
	if check is None:
		qry="INSERT INTO chemical_disease(chemical_index,disease_index) VALUES (%s, %s)"
		try:
			cur.execute(qry,(chemind,disind))
			con.commit()
		except:
			con.rollback()

def update_chem_altid(chemind,altid):
	chemind=str(chemind)
	altid=str(altid)
	qry="UPDATE chemical SET Alternate_id = %s WHERE chemical_index = %s"
	try:
		cur.execute(qry,(altid,chemind))
		con.commit()
	except:
		con.rollback()
def update_chem_altikey(chemind,altikey):
	chemind=str(chemind)
	altikey=str(altikey)
	qry="UPDATE chemical SET Alternate_InChIKey = %s WHERE chemical_index = %s"
	try:
		cur.execute(qry,(altikey,chemind))
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
                        cur.execute(qry,(ikey,altikey,chemname,cid,cas,chembl,altid,smiles))
                        con.commit()
                except:
                        con.rollback()


diseasename_index=get_disease()

S=open('../output/therapeutic_chemical_disease_index.csv',"w")
notfound_file='../output/chemicals_notfound_from_chemdisease.txt'
NF=open(notfound_file,"w")
notfound_file2='../output/chemicals_notfound_from_chemdisease2.txt'
NF2=open(notfound_file2,"w")
for line in open('../therapeuticChemical_disease_CTD.txt',"r").xreadlines():
	line=line.strip().split("\t")
	diseasename_lower=str(line[0])
	chemicalname=str(line[1])
	chemId=str(line[2])
	cas=str(line[3])
	DiseaseName=str(line[4])
	disease_MESH=str(line[5])
	disind=None
	try:
		dis=diseasename_index[diseasename_lower]
		if dis is not None:
			disind=str(dis)
	except:
		None

	chemind=None
	ikey=None
	try:
	#try to find chemical index using CAS, chemicalname, and CTD id
		chem=get_chemical_index_by_CAS(cas)
		if chem is not None:
			chemind=str(chem[0])
		else:
			try:
				ikey=get_InChIKey_by_name(chemicalname)	#pubchem process
				chem=get_chemical_index_by_InChIKey(ikey) #MySQL process
				if chem is not None:
					#alternate InChIKey found from PubChem
					chemind=str(chem[0])
					update_chem_altikey(chemind,ikey) #update alternative InChIKey
				else:
					try:
						cid=get_CID_by_CTD(chemId) #pubchem process
						chem=get_chemical_index_by_CID(cid)
						if chem is not None:
							#alternative CID found in PubChem
							chemind=str(chem[0])
							altid="CID:"+cid
							update_chem_altid(chemind,altid)
						else:
							try:
								altikey=get_InChIKey_by_CAS(cas)
								if altikey is not None:
									#alternative InChIKey found from PubChem
									chem=get_chemical_index_by_InChIKey(ikey) #MySQL process
									if chem is not None:
										#chemical index searchable by ikey in DB
										chemind=str(chem[0])
										update_chem_altikey(chemind,ikey) #update alternative InChIKey
									else:
										#chemical index Not searchable by ikey in DB
										cid=get_CID_by_InChIKey(ikey)
										if cid is not None:
											#alternate CID found in PubChem
											chem=get_chemical_index_by_CID(cid)
											if chem is not None:
												#chemical index searchable by CID in DB
												chemind=str(chem[0])
												altid="CID:"+cid
												update_chem_altid(chemind,altid)
											else:
												NF.write(disind+"\t"+diseasename_lower+"\t"+\
								chemicalname+"\t"+chemId+"\t"+cas+"\t"+disease_MESH+"\t"+altikey+"\t"+cid+"\n")
												
										
							except:
								None
					except:
						None
			except:
				None
		
	except:
		None


	if chemind is not None:
		if disind is not None:
			S.write(chemind+", "+disind+"\n")
			insert_chem_disease(chemind,disind)
		else:
			None
			print "%s, %s\tDisease Not found"%(diseasename_lower,disease_MESH)
	else:
		if disind is not None:
			NF2.write(disind+"\t"+ikey+"\t"+cas+"\t"+chemId)


for line in open(notfound_file, "r").xreadlines():
	line=line.strip().split("\t")
	disind=str(line[0])
	diseasename=str(line[1])
	chemname=str(line[2])
	chemId=str(line[3])
	cas=str(line[4])
	disease_mesh=str(line[5])
	ikey=str(line[6])
	cid=str(line[7])
	smi=None
	if len(ikey)>25:
		cid=get_CID_by_InChIKey(ikey)
		smi=get_canonicalsmiles=by_InChIKey(ikey)
		if len(cas)<3:
			cas=get_CAS_by_InChIKey(ikey)
		
	elif len(cas)>2:
		ikey=get_InChIKey_by_CAS(cas)
		cid=get_CID_by_InChIKey(ikey)
		smi=get_canonicalsmiles_by_InChIKey(ikey)
	
	if len(ikey)>25:
		if smi is not None:
			altikey=None
			if len(chemname)<2:
				chemname=None
			chembl=None
			altid=None
			if len(chemId)>2:
				altid=chemId
			insert_chemical(ikey,altikey,chemname,cid,cas,chembl,altid,smi)
	if disind is not None:
		#additional chemical-disease pairs
		chemind=get_chemical_index_by_InChIKey(ikey)
		S.write(chemind+", "+disind+"\n")
		insert_chem_disease(chemind,disind)

for line in open(notfound_file2, "r").xreadlines():
	line=line.strip().split("\t")
	disind=str(line[0])
	ikey=str(line[1])
	cas=str(line[4])
	chemId=str(line[3])
	cid=None
	smi=None
	chemname=None
	chembl=None
	altid=None
	if len(ikey)>25:
		cid=get_CID_by_InChIKey(ikey)
		smi=get_canonicalsmiles=by_InChIKey(ikey)
		if len(cas)<3:
			cas=get_CAS_by_InChIKey(ikey)
		
	elif len(cas)>2:
		ikey=get_InChIKey_by_CAS(cas)
		cid=get_CID_by_InChIKey(ikey)
		smi=get_canonicalsmiles_by_InChIKey(ikey)
	
	if len(ikey)>25:
		if smi is not None:
			altikey=None
			if len(chemId)>2:
				altid=chemId
			insert_chemical(ikey,altikey,chemname,cid,cas,chembl,altid,smi)
	if disind is not None:
		#additional chemical-disease pairs
		chemind=get_chemical_index_by_InChIKey(ikey)
		S.write(chemind+", "+disind+"\n")
		insert_chem_disease(chemind,disind)
	
con.close()	
