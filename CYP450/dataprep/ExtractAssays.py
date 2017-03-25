#!/usr/bin/python

#from Bio import SeqIO
import MySQLdb as db
import sys

con = db.connect(host='localhost', user='hlim', passwd='secret!', db='chembl_22',unix_socket = '/var/run/mysqld/mysqld.sock');
cur=con.cursor()
def get_assayInfo_by_tid(tid):
    aids=get_aid_by_tid(tid)
    if aids:
        results=[]
        for aid in aids:
            aid=str(aid[0])
            assayType=get_assayType(aid)
            assayTestType=get_assayTestType(aid)
            confScore=get_assay_confidence(aid)
            result=(aid,assayType,assayTestType,confScore)
            results.append(result)
        return results
    else:
        return None
def get_activityInfo_by_aid(aid):
    activityIds=get_activityId_by_aid(aid)
    if activityIds:
        results=[]
        for act in activityIds:
            act=str(act[0])
            molregno=get_molregno_by_activityId(act)
            relation=get_relation_by_activityId(act)
            value=get_value_by_activityId(act)
            unit=get_unit_by_activityId(act)
            actType=get_type_by_activityId(act)
            result=(molregno,relation,value,unit,actType)
            results.append(result)
        return results
    else:
        return None
        
def get_molregno_by_activityId(actId):
    qry="SELECT molregno FROM activities WHERE activity_id=%s"
    cur.execute(qry,(str(actId),))
    res=cur.fetchone()
    if res:
        return str(res[0])
    else:
        return None
def get_relation_by_activityId(actId):
    qry="SELECT standard_relation FROM activities WHERE activity_id=%s"
    cur.execute(qry,(str(actId),))
    res=cur.fetchone()
    if res:
        return str(res[0])
    else:
        return None
def get_value_by_activityId(actId):
    qry="SELECT standard_value FROM activities WHERE activity_id=%s"
    cur.execute(qry,(str(actId),))
    res=cur.fetchone()
    if res:
        return str(res[0])
    else:
        return None
def get_unit_by_activityId(actId):
    qry="SELECT standard_units FROM activities WHERE activity_id=%s"
    cur.execute(qry,(str(actId),))
    res=cur.fetchone()
    if res:
        return str(res[0])
    else:
        return None
def get_type_by_activityId(actId):
    qry="SELECT standard_type FROM activities WHERE activity_id=%s"
    cur.execute(qry,(str(actId),))
    res=cur.fetchone()
    if res:
        return str(res[0])
    else:
        return None
    
def get_activityId_by_aid(aid):
    qry="SELECT activity_id FROM activities WHERE assay_id=%s"
    cur.execute(qry,(str(aid),))
    res=cur.fetchall()
    if res:
        return res
    else:
        return None

def get_aid_by_tid(tid):
    qry="SELECT assay_id FROM assays WHERE tid=%s"
    cur.execute(qry,(str(tid),))
    res=cur.fetchall()
    if res:
        return res
    else:
        return None

def get_assayType(aid):
    qry="SELECT assay_type FROM assays WHERE assay_id=%s"
    cur.execute(qry,(str(aid),))
    res=cur.fetchone()
    if res:
        res=str(res[0])
        if res=="B":
            return "Binding"
        elif res=="A":
            return "ADME"
        elif res=="F":
            return "Functional"
        else:
            "Unknown"
    else:
        return "Unknown"

def get_assayTestType(aid):
    qry="SELECT assay_test_type FROM assays WHERE assay_id=%s"
    cur.execute(qry,(str(aid),))
    res=cur.fetchone()
    if res:
        return str(res[0])
    else:
        return "Unknown"

def get_assay_confidence(aid):
    qry="SELECT confidence_score FROM assays WHERE assay_id=%s"
    cur.execute(qry,(str(aid),))
    res=cur.fetchone()
    if res:
        return str(res[0])
    else:
        return None
def get_chemInfo_by_molregno(molregno):
    ikey=get_ikey_by_molregno(molregno)
    smiles=get_smiles_by_molregno(molregno)
    chiral=get_chirality_by_molregno(molregno)
    prodrug=get_prodrug_by_molregno(molregno)
    result=(ikey,smiles,chiral,prodrug)
    return result

def get_chirality_by_molregno(molregno):
    qry="SELECT chirality FROM molecule_dictionary WHERE molregno=%s"
    cur.execute(qry,(str(molregno),))
    res=cur.fetchone()
    if res:
        res=int(res[0])
        if res==0:
            return "Racemic"
        elif res==1:
            return "Single_Stereoisomer"
        elif res==2:
            return "Achiral"
        else:
            return "Unknown"
    else:
        return "Unknown"
def get_prodrug_by_molregno(molregno):
    qry="SELECT prodrug FROM molecule_dictionary WHERE molregno=%s"
    cur.execute(qry,(str(molregno),))
    res=cur.fetchone()
    if res:
        res=int(res[0])
        if res==0:
            return "No"
        elif res==1:
            return "Yes"
        else:
            return "Unknown"
    else:
        return "Unknown"
def get_ikey_by_molregno(molregno):
    qry="SELECT standard_inchi_key FROM compound_structures WHERE molregno=%s"
    cur.execute(qry,(str(molregno),))
    res=cur.fetchone()
    if res:
        res=str(res[0])
        return res
    else:
        return "Unknown"
def get_smiles_by_molregno(molregno):
    qry="SELECT canonical_smiles FROM compound_structures WHERE molregno=%s"
    cur.execute(qry,(str(molregno),))
    res=cur.fetchone()
    if res:
        res=str(res[0])
        return res
    else:
        return "Unknown"

def count_component_by_tid(tid):
    qry="SELECT tid,component_id FROM target_components WHERE tid=%s"
    cur.execute(qry,(str(tid),))
    res=cur.fetchall()
    return res
def get_protInfo_by_tid(tid):
    compId=get_component_id_by_tid(tid)
    qry="SELECT accession, description, sequence FROM component_sequences WHERE component_id=%s"
    cur.execute(qry,(str(compId),))
    res=cur.fetchone()
    if res:
        accession=str(res[0])
        desc=str(res[1])
        seq=str(res[2])
        result=(tid,accession,desc,seq)
        return result
    else:
        return "Unknown"
def get_component_id_by_tid(tid):
    qry="SELECT component_id FROM target_components WHERE tid=%s"
    cur.execute(qry,(str(tid),))
    res=cur.fetchone()
    if res:
        res=str(res[0])
        return res
    else:
        return "Unknown"
#ainfo=get_assayInfo_by_tid("55")
#print ainfo
#actinfo=get_activityInfo_by_aid("4002")
#print actinfo
#drug=get_chemInfo_by_molregno("115")
#print drug
#sys.exit()  #testing line

assayFile="./CYP450_Assays.csv"
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
AF.write("AssayId, AssayType, TestType, ConfidenceScore, Molregno, TargetId, ActivityType, ValueRelation, ActivityValue, ActivityUnit\n")
for data in dataset:
    line="%s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n" %(data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9])
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
