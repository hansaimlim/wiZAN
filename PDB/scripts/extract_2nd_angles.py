#!/usr/bin/python
import sys
import string
import numpy as np
from Bio.PDB.DSSP import *
import glob
#print(glob.glob("/home/hlim/dssp/2ov*"))
#sys.exit()
outfile='./SeqID30_secondary_str_angle_statistics.txt'
'''
H	Alpha helix (4-12)
B	Isolated beta-bridge residue
E	Strand
G	3-10 helix
I	Pi helix
T	Turn
S	Bend
-	None


Relative accessibility=residue_accessibility / max(residue_accessibility)
'''
h_phi=[]
h_psi=[]
b_phi=[]
b_psi=[]
e_phi=[]
e_psi=[]
g_phi=[]
g_psi=[]
i_phi=[]
i_psi=[]
t_phi=[]
t_psi=[]
s_phi=[]
s_psi=[]
n_phi=[]
n_psi=[]
#########

ciffiles=glob.glob("/home/hlim/PDB/SeqID30/*.cif")
filecount=len(ciffiles)
count=0
for cif in ciffiles:

    f=string.replace(cif,'PDB/SeqID30','dssp')
    f=string.replace(f,'cif','dssp')
    try:
        dssp,keys=make_dssp_dict(f)
    except:
        continue
    keys=keys[1:-1] #Skip the first and last residue as their angles are half-missing
    for key in keys:
        
        ds=dssp[key]
        aa=str(ds[0])
        ss=str(ds[1])
        acc=float(ds[2])
        psi=float(ds[3])
        phi=float(ds[4])
        if (aa=='G') or (aa=='P'):
            #skip glycine or proline
            continue
        if ss=='H':
            h_phi.append(phi)
            h_psi.append(psi)
        elif ss=='B':
            b_phi.append(phi)
            b_psi.append(psi)
        elif ss=='E':
            e_phi.append(phi)
            e_psi.append(psi)
        elif ss=='G':
            g_phi.append(phi)
            g_psi.append(psi)
        elif ss=='I':
            i_phi.append(phi)
            i_psi.append(psi)
        elif ss=='T':
            t_phi.append(phi)
            t_psi.append(psi)
        elif ss=='S':
            s_phi.append(phi)
            s_psi.append(psi)
        else:
            n_phi.append(phi)
            n_psi.append(psi)
    count+=1
    if count % 10000 ==0:
        print("%s files done. Total %s files."%(str(count),str(filecount)))



S=open(outfile,"w")
S.write("StructureType\tmin.Psi\tmax.Psi\tmean.Psi\tstd.Psi\tmin.Phi\tmax.Phi\tmean.Phi\tstd.Phi\tCount\n")

S.write("AlphaHelix\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n"%( np.amin(h_psi),np.amax(h_psi),np.mean(h_psi),np.std(h_psi),np.amin(h_phi),np.amax(h_phi),np.mean(h_phi),np.std(h_phi),str(len(h_psi))))
S.write("BetaBridge\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n"%( np.amin(b_psi),np.amax(b_psi),np.mean(b_psi),np.std(b_psi),np.amin(b_phi),np.amax(b_phi),np.mean(b_phi),np.std(b_phi),str(len(b_psi))))
S.write("Strand\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n"%( np.amin(e_psi),np.amax(e_psi),np.mean(e_psi),np.std(e_psi),np.amin(e_phi),np.amax(e_phi),np.mean(e_phi),np.std(e_phi),str(len(e_psi))))
S.write("3-10Helix\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n"%( np.amin(g_psi),np.amax(g_psi),np.mean(g_psi),np.std(g_psi),np.amin(g_phi),np.amax(g_phi),np.mean(g_phi),np.std(g_phi),str(len(g_psi))))
S.write("PiHelix\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n"%( np.amin(i_psi),np.amax(i_psi),np.mean(i_psi),np.std(i_psi),np.amin(i_phi),np.amax(i_phi),np.mean(i_phi),np.std(i_phi),str(len(i_psi))))
S.write("Turn\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n"%( np.amin(t_psi),np.amax(t_psi),np.mean(t_psi),np.std(t_psi),np.amin(t_phi),np.amax(t_phi),np.mean(t_phi),np.std(t_phi),str(len(t_psi))))
S.write("Bend\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n"%( np.amin(s_psi),np.amax(s_psi),np.mean(s_psi),np.std(s_psi),np.amin(s_phi),np.amax(s_phi),np.mean(s_phi),np.std(s_phi),str(len(s_psi))))
S.write("Others\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n"%( np.amin(n_psi),np.amax(n_psi),np.mean(n_psi),np.std(n_psi),np.amin(n_phi),np.amax(n_phi),np.mean(n_phi),np.std(n_phi),str(len(n_psi))))
