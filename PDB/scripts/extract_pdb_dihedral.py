#!/usr/bin/python

import os
import re
import sys
import math
import subprocess
from multiprocessing import Pool
from multiprocessing import freeze_support, cpu_count
#from Bio.PDB import *
import numpy as np
sys.path.append('/home/hlim/biopython/')
from Bio import PDB as pdb
from Bio import SeqIO as seqio
###################################################### Amino Acid Letter Codes
global amino_acid_3to1
amino_acid_3to1={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLU":"E","GLN":"Q","GLY":"G",\
"HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W",\
"TYR":"Y","VAL":"V","SEC":"U","PYL":"O","GLX":"Z","ASX":"X","UNK":"X","MSE":"X","HCY":"X","HSE":"X",\
"NLE":"X","NVA":"X","ORN":"X","PEN":"X","pGLU":"X","Lys(Dnp)":"X","pTHR":"T","pSER":"S","pTYR":"Y",\
"CIT":"X"}
#global amino_acid_1to3
#for key in amino_acid_3to1:
#    amino_acid_1to3[amino_acid_3to1[key]]=key
###################################################### Amino Acid Letter Codes

def get_cif_files(input_dir):
    filenames=[]
    for fi in os.listdir(input_dir):
        if fi.endswith(".cif"):
            filename=os.path.join(input_dir, fi)
            filenames.append(filename)
    return filenames
def ramachandran_type(residue, next_residue) :
    """Expects Bio.PDB residues, returns ramachandran 'type'

    If this is the last residue in a polypeptide, use None
    for next_residue.

    Return value is a string: "General", "Glycine", "Proline"
    or "Pre-Pro".
    """
    if residue.resname.upper()=="GLY" :
        return "Glycine"
    elif residue.resname.upper()=="PRO" :
        return "Proline"
    elif next_residue is not None \
    and next_residue.resname.upper()=="PRO" :
        #exlcudes those that are Pro or Gly
        return "Pre-Pro"
    else :
        return "General"    
def rad_to_degree(rad):
    if rad is None:
        return None
    else:
        angle=rad*float(180.0)/np.pi
        if angle>180.0:
            angle=angle-float(360.0)
        elif angle<-180.0:
            angle=angle+float(360.0)
        return angle

def get_cos_angle(norm1,norm2):
    return (float(180.0)/(np.pi))*np.arccos(np.dot(norm1,norm2)/(np.linalg.norm(norm1)*np.linalg.norm(norm2)))

def extract_Ca(cif_file,Dh_path):    
    namepiece=re.search( r'([a-z0-9.]+).cif', cif_file, re.M)
    Dihedral_file=Dh_path+namepiece.group(1)+"_dihedral.tsv"
    print("Processing file: "+cif_file)
    #Information collected from the CIF file
    parser = pdb.MMCIFParser()
    structure = parser.get_structure(namepiece.group(1).upper, cif_file)
    polypeptides = pdb.CaPPBuilder().build_peptides(structure)
    print("Calculating dihedrarl angles for %s ..." %cif_file)
    DH=open(Dihedral_file,"w")
    DH.write("Chain:Residue\tNum\tPhi\tPsi\tRamachandran_Type\n")
    for model in structure :
        for chain in model :
            print("Chain %s" % str(chain.id))
            polypeptides = pdb.CaPPBuilder().build_peptides(chain)
            for poly_index, poly in enumerate(polypeptides) :
                phi_psi = poly.get_phi_psi_list()
                seq=poly.get_sequence()
                print(seq)
             #   seqio.write(seq,'./outputseq.txt','fasta')
                for res_index, residue in enumerate(poly) :
                     phi, psi = phi_psi[res_index]
                     if phi and psi :
                         #Don't write output when missing an angle
                        DH.write("%s:%s\t%i\t%f\t%f\t%s\n" \
                         % (str(chain.id), residue.resname, residue.id[1], rad_to_degree(phi), rad_to_degree(psi), ramachandran_type(residue, poly[res_index+1])))
    DH.close()

if __name__ == '__main__':
    Args=sys.argv[1:]
    if len(Args)<2:
        print("Usage: python3 extract_pdb.py input_dir dihedral_out_dir")
        print("Example: python3 extract_pdb.py /home/cif_dir/ /home/dihedral_dir/")
        sys.exit()
    filenames=get_cif_files(Args[0])
#    print(filenames)
#    print(filenames[0])
#    print(filenames[1])
#    print(filenames[2])
#    namepiece=re.search( r'([a-z0-9]+).cif', filenames[0], re.M)
#    print(namepiece.group(1))
#    sys.exit()
#    extract_Ca(filenames[0],Args[1],Args[2])
    inputs=[]
    for filename in filenames:
        inputs.append((filename,Args[1]))
#    print(filenames.index("./SeqID30/1usc.cif"))
#    print(filenames.index("./SeqID30/3k3r.cif"))
#    print(filenames.index("./SeqID30/3r4k.cif"))
#    print(filenames[23267])
#    norm1=np.array([-0.4765,0.8713,1.3115],dtype=np.float32)
#    norm2=np.array([0.8531,-1.9017,-0.2123],dtype=np.float32)
#    print(norm1-norm2)
#    print(get_cos_angle(norm1,norm2))
#    sys.exit()
    freeze_support()
    with Pool(cpu_count()-2) as pool:
        try:
            R=pool.starmap(extract_Ca,inputs)
        except BaseException as e:
            print("Error occurred:\n"+str(e)+"\n")
    
