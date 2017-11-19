#!/usr/bin/python
import sys
import string
import numpy as np
import glob
from multiprocessing import Pool
from multiprocessing import freeze_support, cpu_count
#from Bio.PDB import *
sys.path.append('/home/hlim/biopython/')
from Bio.PDB.DSSP import *
from Bio import PDB as pdb
from Bio import SeqIO as seqio
from Bio.PDB.Polypeptide import three_to_one, one_to_three
#print(glob.glob("/home/hlim/dssp/2ov*"))
#sys.exit()

def accessibility_class(residue, accessibility):
    #get solvent accessibility class
    #use relative accessibility.
    #acc>=0.95 (2), 0.95>acc>=0.05 (1), 0.05>acc>0 (0)
    Type='Miller' #Miller or Wilke type available
    resmax=residue_max_acc[Type]
    try:
        rel_acc=float(accessibility)/float(resmax[one_to_three(residue)])
    except:
        return ("NA","NA")
#    print(rel_acc)
    if rel_acc>=0.95:
        return (rel_acc,2)
    elif rel_acc>=0.05:
        return (rel_acc,1)
    else:
        return (rel_acc,0)

def angle_class(phi,psi):
    #get dihedral angle class
    #find the closest secondary structure by mean of phi, psi angles
    #does not consider other factors than means
    phi=float(phi)
    psi=float(psi)
    Type='withoutGP' #withGP available. angle data with/without glycine and proline; alldssp for total 130884 structures
    alldssp =    {'H':{'phi':-64.65,'psi':-39.75},'B':{'phi':-97.76,'psi':122.85},'E':{'phi':-110.89,'psi':122.23},'G':{'phi':-65.90,'psi':15.63},\
                  'I':{'phi':-78.71,'psi':-42.23},'T':{'phi':-42.40,'psi':4.50},'S':{'phi':-68.02,'psi':38.36},'-':{'phi':-71.43,'psi':103.70}}
    withGP =     {'H':{'phi':-64.49,'psi':-39.935},'B':{'phi':-96.792,'psi':123.607},'E':{'phi':-110.914,'psi':122.162},'G':{'phi':-66.901,'psi':-15.527},\
                  'I':{'phi':-79.19,'psi':-42.151},'T':{'phi':-42.09,'psi':4.978},'S':{'phi':-68.322,'psi':39.994},'-':{'phi':-69.98,'psi':104.746}}
    withoutGP =  {'H':{'phi':-64.809,'psi':-39.941},'B':{'phi':-103.507,'psi':129.329},'E':{'phi':-115.192,'psi':126.233},'G':{'phi':-70.503,'psi':-13.964},\
                  'I':{'phi':-80.252,'psi':-42.37 },'T':{'phi':-62.26  ,'psi':5.432  },'S':{'phi':-86.235 ,'psi':49.454 },'-':{'phi':-78.27 ,'psi':108.138 }}

    if Type=='withoutGP':
        angles=withoutGP
    elif Type=='alldssp':
        angles=alldssp
    elif Type=='withGP':
        angles=withGP
    else:
        print("Error: choose correct angle type.")
        sys.exit()

    distances={}
    for ss,angle in angles.items():
        xdist1=np.absolute(phi-float(angle['phi']))
        ydist1=np.absolute(psi-float(angle['psi']))
        if phi>angle['phi']:
            xdist2=np.absolute(phi-float(360)-float(angle['phi']))
        else:
            xdist2=np.absolute(phi+float(360)-float(angle['phi']))
        if psi>angle['psi']:
            ydist2=np.absolute(psi-float(360)-float(angle['psi']))
        else:
            ydist2=np.absolute(psi+float(360)-float(angle['psi']))

        phidist=np.amin([xdist1,xdist2])
        psidist=np.amin([ydist1,ydist2])
        distance=np.sqrt(np.square(phidist)+np.square(psidist))
        distances[ss]=distance
#    print(distances)
    distlist=sorted(distances.items(),key=lambda x: x[1])
    return distlist[0][0]
    
        

if __name__ == '__main__':
#    Args=sys.argv[1:]
#    angle_class(0,0)
#    angle_class(180,-180)
#    accessibility_class('A',90)

    ciffiles=glob.glob("/home/hlim/PDB/SeqID30/*.cif")
    filecount=len(ciffiles)
    count=0
    for cif in ciffiles:
        sequence=''
        f=string.replace(cif,'PDB/SeqID30','dssp')
        f=string.replace(f,'cif','dssp')
        outfile=string.replace(cif,'SeqID30','SeqID30_secondary')
        outfile=string.replace(outfile,'cif','tsv')
        pdbid=string.replace(cif,'/home/hlim/PDB/SeqID30/','')
        pdbid=string.replace(pdbid,'.cif','')
        OUT=open(outfile,"w")
        
        try:
            dssp,keys=make_dssp_dict(f)
        except:
            OUT.close()
            continue
        OUT.write("AA\t2ndStr\tRel_acc\tAcc_class\tPhi\tPsi\tAngle_class_closest\n")
        keys=keys[1:-1] #Skip the first and last residue as their angles are half-missing
        for key in keys:

            ds=dssp[key]
            aa=str(ds[0])
            sequence=sequence+aa
            ss=str(ds[1])
            acc=float(ds[2])
            phi=float(ds[3])
            psi=float(ds[4])
            accClass=accessibility_class(aa,acc)
            angleClass=angle_class(phi,psi)
            try:
                rel_acc=float(accClass[0])
                OUT.write("%s\t%s\t%.4f\t%s\t%s\t%s\t%s\n"%(aa,ss,rel_acc,accClass[1],str(phi),str(psi),str(angleClass)))
            except:
                OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(aa,ss,accClass[0],accClass[1],str(phi),str(psi),str(angleClass)))
            
        OUT.close()
        outfile_sequence=string.replace(outfile,'tsv','fas')
        OUTSEQ=open(outfile_sequence,"w")
        OUTSEQ.write(">%s Sequence from DSSP module\n"%pdbid)
        OUTSEQ.write(sequence+"\n")
        count+=1
        if count % 10000 ==0:
            print("%s files done. Total %s files."%(str(count),str(filecount)))


