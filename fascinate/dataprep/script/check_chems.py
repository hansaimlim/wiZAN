#!/usr/bin/python

ikeyfile='../uniq_ikeys.txt'
infofile='../output/chem_info_from_PubChem.tsv'
cmapfile='../cMap_InChIKey_drugname.tsv'
ZCDfile='../../../ZINC_ChEMBL_DrugBank/chem_chem/ZINC_ChEMBL_DrugBank_chemicals.tsv'
FDAfile='../../../ZINC_ChEMBL_DrugBank/chem_chem/FDA_InChIKey_compoundname.tsv'
ikeys=[]
ikeys_found=[]
ikey_name={}
ikey_smiles={}
for line in open(ikeyfile, "r").xreadlines():
	line=line.strip().split("\t")
	ikey=str(line[0])
	ikeys.append(ikey)
for line in open(cmapfile, "r").xreadlines():
	line=line.strip().split("\t")
	ikey=str(line[0])
	drugname=str(line[1])
	try:
		ikey_name[ikey]
	except:
		ikey_name[ikey]=drugname
for line in open(ZCDfile,"r").xreadlines():
	line=line.strip().split("\t")
	ind=str(line[0])
	ikey=str(line[1])
	smiles=str(line[2])
	try:
		ikey_smiles[ikey]
	except:
		ikey_smiles[ikey]=smiles
for line in open(FDAfile, "r").xreadlines():
	line=line.strip().split("\t")
	ikey=str(line[0])
	drugname=str(line[1])
	try:
		ikey_name[ikey]
	except:
		ikey_name[ikey]=drugname

for line in open(infofile, "r").xreadlines():
	line=line.strip().split("\t")
	ikey=str(line[0])
	cid=str(line[1])
	synonym=str(line[2])
	cas=str(line[3])
	smiles=str(line[4])
	ikeys_found.append(ikey)

notfile='../output/ikey_infoNotFound.txt'
addfile1='../output/ikey_infoFound_drugname.txt'
addfile2='../output/ikey_infoFound_smilesOnly.txt'
NF=open(notfile,"w")
A1=open(addfile1,"w")
A2=open(addfile2,"w")
for ikey in ikeys:
	if (ikey in ikeys_found):
		continue
	else:
		#the given inchikey not found from pubchem
		try:
			drugname=ikey_name[ikey]
			smiles=ikey_smiles[ikey]
			line="%s\t%s\t%s"%(ikey,drugname,smiles)
			A1.write(line+"\n")
		except:
			try:
				smiles=ikey_smiles[ikey]
				line="%s\t%s"%(ikey,smiles)
				A2.write(line+"\n")
			except:
				NF.write(ikey+" was not found\n")

