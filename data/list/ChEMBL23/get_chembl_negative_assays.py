#!/usr/bin/python
import MySQLdb as db
import sys

def get_negative_activities():
    #collect activity information from ChEMBL db
    qry="SELECT act.assay_id, act.molregno, ay.tid, act.standard_value, act.standard_units\
         FROM chembl_23.activities act INNER JOIN chembl_23.assays ay ON act.assay_id=ay.assay_id\
         WHERE act.potential_duplicate = 0 AND act.standard_type = 'IC50' AND act.standard_value > 20000\
         AND ay.confidence_score=9"
    cur.execute(qry)
    result=cur.fetchall()
    return result
def split_seq(seq, num):
    return [ seq[start:start+num] for start in range(0, len(seq), num) ]
def get_compound(molregno):
    val=(str(molregno),)
    qry="SELECT cstr.standard_inchi_key, cstr.canonical_smiles FROM chembl_23.compound_structures cstr\
         WHERE cstr.molregno=%s"
    cur.execute(qry,val)
    result=cur.fetchone()
    return result
def get_target(tid):
    val=(str(tid),)
    qry="SELECT cseq.accession, cseq.sequence, cseq.description, csyn.component_synonym FROM component_sequences cseq\
         INNER JOIN target_components tcomp ON cseq.component_id=tcomp.component_id\
         INNER JOIN component_synonyms csyn ON tcomp.component_id=csyn.component_id\
         WHERE tcomp.tid=%s AND cseq.component_type='PROTEIN'"
    cur.execute(qry,val)
    result=cur.fetchone()
    return result
def get_cheminfo(molregno):
    val=(str(molregno),)
    qry="SELECT chembl_id, pref_name, max_phase, therapeutic_flag, availability_type from molecule_dictionary\
         WHERE molregno=%s"
    cur.execute(qry,val)
    result=cur.fetchone()
    return result



passwd=raw_input("MySQL password?")
con=db.connect(host='127.0.0.1',user='hlim',passwd=passwd,db='chembl_23');
cur=con.cursor()

chem_prot_file='./ChEMBL23_negative_assays_20uM.txt'
chem_prot_index='ChEMBL23_negative_assays_index_20uM.txt'
#chem_prot_file='./ChEMBL23_chem_prot_active_IC50_10uM.tsv'
#chem_prot_index='./ChEMBL23_chem_prot.csv'
chemInfo_file='./ChEMBL23_chemInfo.tsv'
protInfo_file='./ChEMBL23_protInfo.tsv'
prot_fasta_file='./ChEMBL23_prot_seq.fas'
print "Start parsing ChEMBL23 database..."
activities=get_negative_activities()
molregnos=[]
tids=[]
for activity in activities:
    aid=activity[0]
    molregno=str(activity[1])
    tid=str(activity[2])
    ic50=activity[3]
    molregnos.append(molregno)
    tids.append(tid)
    

molregnos=list(set(molregnos))
tids=list(set(tids))

idx2ikey={}
ikey2idx={}
mol2idx={}
idx2acc={}
acc2idx={}
tid2idx={}
protinfo=[]
idx2chembl={}
chembl2idx={}

mol2idx={}
idx2mol={}
maxchemIdx=1
with open(chemInfo_file,"r") as inf:
    next(inf)
    for line in inf:
        line=line.strip().split("\t")
        chemindex=int(line[0])
        molregno=str(line[1])
        ikey=str(line[2])
        smiles=str(line[3])
        mol2idx[molregno]=chemindex
        idx2mol[chemindex]=molregno
        idx2ikey[chemindex]=ikey
        maxchemIdx=chemindex
with open(protInfo_file,"r") as inf:
    next(inf)
    for line in inf:
        line=line.strip().split("\t")
        protindex=int(line[0])
        tid=str(line[1])
        acc=str(line[2])
        syn=str(line[3])
        desc=str(line[4])
        idx2acc[protindex]=acc
        acc2idx[acc]=protindex
        tid2idx[tid]=protindex
with open('./ChEMBL23_chemInfo_detailed.tsv',"a") as chemdetail:
    for mol in molregnos:
        mol=str(mol)
        try:
            idx=mol2idx[mol]
        except:
            #new chemical
            try:
                comprec=get_compound(mol) #(ikey,smiles)
                chemrec=get_cheminfo(mol) #(chemblID,prefname,maxphase,therapeutic,availability)
            except:
                continue
            try:
                ikey=str(comprec[0])
                smi=str(comprec[1])
                chembl=str(chemrec[0])
            except:
                continue
            maxchemIdx+=1 #new chemical index
            mol2idx[mol]=maxchemIdx
            idx2mol[maxchemIdx]=mol
            try:
                pref_name=str(chemrec[1])
            except:
                pref_name="None"
            try:
                max_phase=str(chemrec[2])
            except:
                max_phase="None"
            therapeutic_flag=str(chemrec[3])
            try:
                availability_type=str(chemrec[4])
            except:
                availability_type="None"
            chemdetail.write(str(maxchemIdx)+"\t"+chembl+"\t"+str(mol)+"\t"+ikey+"\t"+pref_name+"\t"+max_phase+"\t"+therapeutic_flag+"\t"+availability_type+"\t"+smi+"\n")
       

 

with open(chem_prot_index,"w") as fid:
    for activity in activities:
        aid=activity[0]
        molregno=str(activity[1])
        tid=str(activity[2])
        ic50=activity[3]
        try:
            chemIdx=mol2idx[molregno]
            protIdx=tid2idx[tid]
        except:
            continue
        fid.write("%s, %s, -1\n"%(str(chemIdx),str(protIdx)))
        

with open(chem_prot_file,"w") as fid:
    fid.write("ChemIdx\tProtIdx\tIC50(nM)\tMolregno\tInChIKey\tTID\tAccession\n")
    for activity in activities:
        aid=activity[0]
        molregno=str(activity[1])
        tid=str(activity[2])
        ic50=activity[3]
        try:
            chemIdx=mol2idx[molregno]
            protIdx=tid2idx[tid]
            ikey=idx2ikey[chemIdx]
            acc=idx2acc[protIdx]
        except:
            continue
        fid.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(str(chemIdx),str(protIdx),str(ic50),molregno,str(ikey),tid,str(acc)))
