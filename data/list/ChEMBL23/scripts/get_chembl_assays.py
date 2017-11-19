#!/usr/bin/python
import MySQLdb as db
import sys

def get_activities():
    #collect activity information from ChEMBL db
    qry="SELECT act.assay_id, act.molregno, ay.tid, act.standard_value, act.standard_units\
         FROM chembl_23.activities act INNER JOIN chembl_23.assays ay ON act.assay_id=ay.assay_id\
         WHERE act.potential_duplicate = 0 AND act.standard_type = 'IC50' AND act.standard_value <= 10000\
         AND ay.confidence_score>=8"
    cur.execute(qry)
    result=cur.fetchall()
    return result
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
def split_seq(seq, num):
    return [ seq[start:start+num] for start in range(0, len(seq), num) ]

passwd=raw_input("MySQL password?")
con=db.connect(host='127.0.0.1',user='hlim',passwd=passwd,db='chembl_23');
cur=con.cursor()

chem_prot_file='./temp'
chem_prot_index='./temp'
chemInfo_file='./temp'
protInfo_file='./temp'
prot_fasta_file='./ChEMBL23_prot_seq.fas'
#chem_prot_file='./ChEMBL23_chem_prot_active_IC50_10uM.tsv'
#chem_prot_index='./ChEMBL23_chem_prot.csv'
#chemInfo_file='./ChEMBL23_chemInfo.tsv'
#protInfo_file='./ChEMBL23_protInfo.tsv'
#prot_fasta_file='./ChEMBL23_prot_seq.fas'
print "Start parsing ChEMBL23 database..."
activities=get_activities()
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

chemIdx=1
with open(chemInfo_file,"w") as fid:
#    fid.write("Index\tMolregno\tInChIKey\tCanonicalSMILES\n")
    for mol in molregnos:
        record=get_compound(mol)
        try:
            ikey=record[0]
            smi=record[1]
        except:
            continue
        mol2idx[mol]=chemIdx
        idx2ikey[chemIdx]=ikey
        ikey2idx[ikey]=chemIdx
#        fid.write("%s\t%s\t%s\t%s\n" %(str(chemIdx),str(mol),str(ikey),str(smi)))
        chemIdx+=1
protIdx=1
with open(prot_fasta_file, "a") as fid:
   # fid.write("Index\tTID\tAccession\tSynonym\tDescription\n")
    for tid in tids:
        record=get_target(tid)
        try:
            acc=str(record[0])
            seq=str(record[1])
        except:
            continue
        try:
            desc=str(record[2])
            syn=str(record[3])
        except:
            desc="No Description"
            syn="No Synonym"
        idx2acc[protIdx]=acc
        acc2idx[acc]=protIdx
        tid2idx[tid]=protIdx
        prot=(protIdx,acc,syn,desc,seq)
        fid.write(">%s |%s |%s |%s\n"%(str(protIdx),str(acc),str(syn),str(desc)))
        fasta=split_seq(seq,70)
        for ff in fasta:
            fid.write(ff+"\n")
        protinfo.append(prot)
  #      fid.write("%s\t%s\t%s\t%s\t%s\n"%(str(protIdx),str(tid),acc,syn,desc))
        protIdx+=1
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
#        fid.write("%s, %s, 1\n"%(str(chemIdx),str(protIdx)))
        

with open(chem_prot_file,"w") as fid:
#    fid.write("ChemIdx\tProtIdx\tIC50(nM)\tMolregno\tInChIKey\tTID\tAccession\n")
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
#        fid.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(str(chemIdx),str(protIdx),str(ic50),molregno,str(ikey),tid,str(acc)))
