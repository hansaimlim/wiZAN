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
        angle=float(180.0)/np.pi
        if angle>180.0:
            angle=angle-float(360.0)
        elif angle<-180.0:
            angle=angle+float(360.0)
        return angle

def get_cos_angle(norm1,norm2):
    return (float(180.0)/(np.pi))*np.arccos(np.dot(norm1,norm2)/(np.linalg.norm(norm1)*np.linalg.norm(norm2)))

def get_dihedral(Coord_CA, Coord_C, Coord_N):
    Phi={}
    Psi={}
    for i in Coord_CA:
        i=int(i)
        phi=float(0.0)
        psi=float(0.0)
        N_curr=np.array(Coord_N[i],dtype=np.float32)
        CA=np.array(Coord_CA[i],dtype=np.float32)
        C_curr=np.array(Coord_C[i],dtype=np.float32)
        alphanorm=np.cross((C_curr-CA),(N_curr-CA))
        try:
            C_prev=np.array(Coord_C[i-1],dtype=np.float32)
            phinorm=np.cross((N_curr-CA),(C_prev-N_curr))
            phi=get_cos_angle(alphanorm,phinorm)
            if phi>float(180.0):
                phi=phi-float(360.0)
        except:
            phi='NA'
        try:
            N_next=np.array(Coord_N[i+1],dtype=np.float32)
            psinorm=np.cross((C_curr-CA),(N_next-C_curr))
            psi=get_cos_angle(psinorm,alphanorm)
            if psi>float(180.0):
                psi=psi-float(360.0)
        except:
            psi='NA'
        Phi[i]=phi
        Psi[i]=psi
    return Psi, Phi
def extract_Ca(cif_file,fasta_path,Ca_path,Dh_path):    
    namepiece=re.search( r'([a-z0-9.]+).cif', cif_file, re.M)
    fasta_file=fasta_path+namepiece.group(1)+".fas"
    Ca_file=Ca_path+namepiece.group(1)+"_Ca.tsv"
    Ca_distmat_file=Ca_path+namepiece.group(1)+"_CaMat.tsv"
    Dihedral_file=Dh_path+namepiece.group(1)+"_dihedral.tsv"
    print("Processing file: "+cif_file)
    #Information collected from the CIF file
    Coord_CA={}
    Coord_C={} #Coord_C[amino_acid_index]=[x,y,z]
    Coord_N={}
    amino_acid_list={} #one-letter code amino acid for the position
    resolution=0.0
    R_work=0.0
    R_free=0.0
    seq_length=0
    #Information collected from the CIF file
    
    Switches={} #turn on when needed
    Switches['xyz']=0 #define swtiches
    CA_probability=float(0.0) #should add up to 1 for each alpha-carbon
    with open(cif_file, "r") as inf:
        X_temp,Y_temp,Z_temp,prob_temp=[],[],[],[]
        for line in inf:
            line=line.strip().split()
            if len(line)<1: #empty lines
                continue
            if line[0]=="#":
                #tag closed
                for key in Switches:
                    Switches[key]=0
                continue
            elif line[0]=="_refine.ls_d_res_high":
                #resolution in Angstrom
                try:
                    resolution=float(line[1])
                except:
                    resolution="NA"
                continue
            elif line[0]=="_refine.ls_R_factor_R_work":
                #R-work
                try:
                    R_work=float(line[1])
                except:
                    R_work="NA"
                continue
            elif line[0]=="_refine.ls_R_factor_R_free":
                #R-free
                try:
                    R_free=float(line[1])
                except:
                    R_free="NA"
                continue
            elif line[0]=="_atom_site.pdbx_PDB_model_num":
                #coordinate information swtich on
                Switches['xyz']=1
                continue
            if Switches['xyz']:
                #get coordinate information
                aa_type3=line[5] #3-letter code for amino acids in CAPS
                try:
                    aa_type1=amino_acid_3to1[aa_type3]
                except:
                    continue
                if line[6]!='A':
                    continue
                atom_type=line[3] #CA for alpha carbon; OG for gamma oxygen, etc.
                aa_num=int(line[8])
                x,y,z,prob=line[10:14]
                if (line[3]=='CA') or (line[3]=='C') or (line[3]=='N'):
                    #alpha carbon in protein chain
                    if float(prob)<1.0:
                        #multiple CA probability and coordinates
                        X_temp.append(x)
                        Y_temp.append(y)
                        Z_temp.append(z)
                        prob_temp.append(prob)
                        if np.sum(np.array(prob_temp,dtype=np.float32))==1.0: #finally all probs found
                            x_avg,y_avg,z_avg=0.0,0.0,0.0
                            for i in range(len(X_temp)):
                                x_avg+=float(X_temp[i])*float(prob_temp[i])
                                y_avg+=float(Y_temp[i])*float(prob_temp[i])
                                z_avg+=float(Z_temp[i])*float(prob_temp[i])
                            coord=np.array([x_avg,y_avg,z_avg],dtype=np.float32)
                            X_temp,Y_temp,Z_temp,prob_temp=[],[],[],[]
                            if line[3]=='CA':
                                Coord_CA[aa_num]=coord
                                amino_acid_list[aa_num]=aa_type1
                            elif line[3]=='C':
                                Coord_C[aa_num]=coord
                            elif line[3]=='N':
                                Coord_N[aa_num]=coord
                        else:
                            continue
                    else: #Single probability and coordinates
                        coord=np.array([x,y,z],dtype=np.float32)
                        if line[3]=='CA':
                            Coord_CA[aa_num]=coord
                            aa_type1=amino_acid_3to1[aa_type3]
                            amino_acid_list[aa_num]=aa_type1
                        #clear temporary coordinates
                        elif line[3]=='C':
                            Coord_C[aa_num]=coord
                        elif line[3]=='N':
                            Coord_N[aa_num]=coord
                        X_temp,Y_temp,Z_temp,prob_temp=[],[],[],[]
                          
                else:
                    continue  

    sequence=''
    CA=open(Ca_file,"w+")
    print("Output alpha-carbon trace to: "+Ca_file)
    CA.write("Num\tAA\tX\tY\tZ\n")
    for aa_num in amino_acid_list:
        aa=amino_acid_list[aa_num]
        CA_xyz=Coord_CA[aa_num]
        CA.write(str(aa_num)+"\t"+str(aa)+"\t"+str(CA_xyz[0])+"\t"+str(CA_xyz[1])+"\t"+str(CA_xyz[2])+"\n")
        sequence=sequence+aa
    CA.close()

    seq_length=len(sequence)
    FA_header=">"+namepiece.group(1)+"|len:"+str(seq_length)+"|res:"+str(resolution)+"|Rfree:"+str(R_free)+"|Rwork:"+str(R_work)+"\n"
    FA=open(fasta_file,"w+")
    print("Output FASTA to: "+fasta_file)
    FA.write(FA_header)
    seq_chunks=(sequence[0+i:70+i] for i in range(0, len(sequence), 70))
    for seq in seq_chunks:
        FA.write(seq+"\n")
    FA.close()
    print("Output alpha-carbon distance matrix to: "+Ca_distmat_file)
    CAMAT=open(Ca_distmat_file,"w")
    for i in amino_acid_list: #header
        CAMAT.write("\t"+str(amino_acid_list[i]))
    CAMAT.write("\n")

    for i in amino_acid_list:
        xyz1=np.array(Coord_CA[i],dtype=np.float32)
        CAMAT.write(str(amino_acid_list[i]))
        for j in amino_acid_list:
            xyz2=np.array(Coord_CA[j],dtype=np.float32)
            dist=np.linalg.norm(xyz1-xyz2)
            CAMAT.write("\t"+str(dist))
        CAMAT.write("\n")
    CAMAT.close()
    print("Extracting %s complete"%cif_file)

    parser = MMCIFParser()
    structure = parser.get_structure(namepiece.group(1), cif_file)
    print("Calculating dihedrarl angles for %s ..."%cif_file)
    DH=open(Dihedral_file,"w")
    DH.write("Chain:Residue:ResID\tPhi\tPsi\tRamachandran_Type\n")
    for model in structure :
        for chain in model :
            print "Chain %s" % str(chain.id)
            polypeptides = Bio.PDB.CaPPBuilder().build_peptides(chain)
            for poly_index, poly in enumerate(polypeptides) :
                phi_psi = poly.get_phi_psi_list()
                for res_index, residue in enumerate(poly) :
                     phi, psi = phi_psi[res_index]
                     if phi and psi :
                         #Don't write output when missing an angle
                        DH.write("Chain%s:%s%i\t%f\t%f\t%s\n" \
                         % (str(chain.id), residue.resname, residue.id[1], degrees(phi), degrees(psi), ramachandran_type(residue, poly[res_index+1])))
    DH.close()

if __name__ == '__main__':
    Args=sys.argv[1:]
    if len(Args)<4:
        print("Usage: python3 extract_pdb.py input_dir fasta_out_dir C_alpha_out_dir dihedral_out_dir")
        print("Example: python3 extract_pdb.py /home/cif_dir/ /home/fasta_dir/ /home/CA_dir/ /home/dihedral_dir/")
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
        inputs.append((filename,Args[1],Args[2],Args[3]))
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
    with Pool(cpu_count()-3) as pool:
        try:
            R=pool.starmap(extract_Ca,inputs)
        except BaseException as e:
            print("Error occurred:\n"+str(e)+"\n")
    
