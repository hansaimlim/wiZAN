#!/usr/bin/python

#from Bio import SeqIO
import MySQLdb as db
import sys

con = db.connect(host='localhost', user='hlim', passwd='w31c0m3', db='chembl_22',unix_socket = '/var/run/mysqld/mysqld.sock');
cur=con.cursor()
def get_mdl_by_molregno(molregno):
    molregno=str(molregno)
    qry="SELECT molfile FROM compound_structures WHERE molregno=%s"
    cur.execute(qry,(molregno,))
    res=cur.fetchone()
    if res:
        return res[0]
    else:
        return None
chemfile="./CYP450_chemicals.tsv"
mdlfile="./CYP450_chemicals_MDL.mol"
M=open(mdlfile,"w")
count=0
with open(chemfile,"r") as inf:
    next(inf)
    for line in inf:
        line=line.strip().split("\t")
        molregno=str(line[0])
        mdl=get_mdl_by_molregno(molregno)
        if mdl:
            M.write(str(molregno))
            M.write(str(mdl)+"\n")
        else:
            print molregno
            count+=1
print count
sys.exit()
#ainfo=get_assayInfo_by_tid("55")
#print ainfo
#actinfo=get_activityInfo_by_aid("4002")
#print actinfo
#drug=get_chemInfo_by_molregno("115")
#print drug
#sys.exit()  #testing line

assayFile="./CYP450_Assays.tsv"
targetInfo="./CYP450_Targets.tsv" #comma may be part of description
chemInfo="./CYP450_chemicals.tsv" #comma may be part of SMILES
#domainInfo="./CYP450_domains.csv"
AF=open(assayFile,"w")
TF=open(targetInfo,"w")
CF=open(chemInfo,"w")
#DF=open(domainInfo,"w")
molregnos=[]
domainids=[]

qry="SELECT DISTINCT(tid) FROM target_dictionary WHERE pref_name LIKE 'Cytochrome P450%' AND tax_id=9606 AND target_type='Single Protein'"
cur.execute(qry)
result = cur.fetchall()
human_cyp_tids=[]
for s in result:
    human_cyp_tids.append(str(s[0]))
#print len(human_cyp_tids)
TF.write("TargetID\tAccession\tDescription\tSequence\n")
for tid in human_cyp_tids:
    tid=str(tid)
    prot=get_protInfo_by_tid(tid)
    tid=prot[0]
    accession=prot[1]
    description=prot[2]
    sequence=prot[3]
    TF.write(tid+"\t"+accession+"\t"+description+"\t"+sequence+"\n")
    
dataset = [] #(aid, type, invitro, conf_score, compound_id, target_id, standard_type, standard_relation, standard_value, standard_units)
for tid in human_cyp_tids:
    assayInfo=get_assayInfo_by_tid(tid)
    for ai in assayInfo:
        aid = str(ai[0])
        activityInfo = get_activityInfo_by_aid(aid)
        for s in activityInfo:
            molregnos.append(s[0])
            dataset.append((ai[0], ai[1], ai[2], ai[3], s[0], tid, s[4], s[1], s[2], s[3]) )
AF.write("AssayId\tAssayType\tTestType\tConfidenceScore\tMolregno\tTargetId\tActivityType\tValueRelation\tActivityValue\tActivityUnit\n")
for data in dataset:
    line="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9])
    AF.write(line)
#print len(dataset)

molregnos=list(set(molregnos))
CF.write("Molregno\tStandardInChIKey\tChirality\tIsProdrug\tCanonicalSMILES\n")
for molregno in molregnos:
    chemInfo=get_chemInfo_by_molregno(str(molregno))
    ikey=str(chemInfo[0])
    smiles=str(chemInfo[1])
    chiral=str(chemInfo[2])
    prodrug=str(chemInfo[3])
    CF.write(str(molregno)+"\t"+ikey+"\t"+chiral+"\t"+prodrug+"\t"+smiles+"\n")
