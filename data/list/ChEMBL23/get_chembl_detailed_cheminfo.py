#!/usr/bin/python
import MySQLdb as db
import sys

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

chemInfo_file='./ChEMBL23_chemInfo.tsv'
protInfo_file='./ChEMBL23_protInfo.tsv'
new_cheminfo='./ChEMBL23_chemInfo_detailed.tsv'
print "Start parsing ChEMBL23 database..."

S=open(new_cheminfo,"w")
S.write("ChemIndex\tChEMBL_id\tMolregno\tInChIKey\tPreferred_name\tMax_phase\tIsTherapeutic\tAvailability\tCanonical_SMILES\n")
with open(chemInfo_file,"r") as inf:
    next(inf)
    for line in inf:
        line=line.strip().split("\t")
        chemidx=str(line[0])
        molregno=str(line[1])
        ikey=str(line[2])
        smi=str(line[3])
        record=get_cheminfo(molregno)
        chembl_id=str(record[0])
        try:
            pref_name=str(record[1])
        except:
            pref_name="None"
        try:
            max_phase=str(record[2])
        except:
            max_phase="None"
        therapeutic_flag=str(record[3])
        try:
            availability_type=str(record[4])
        except:
            availability_type="None"
        S.write(chemidx+"\t"+chembl_id+"\t"+molregno+"\t"+ikey+"\t"+pref_name+"\t"+max_phase+"\t"+therapeutic_flag+"\t"+availability_type+"\t"+smi+"\n")
