#!/usr/bin/python

import sys
import numpy as np
import pandas as pd
from theano import config


extfile='./cyp450_chems_ExtFP.csv'
pubfile='./cyp450_chems_PubChemFP.csv'
subfile='./cyp450_chems_SubstructFP.csv'
#protIdx={'CYP1A2':1,'CYP2C9':2,'CYP2C19':3,'CYP2D6':4,'CYP3A4':5} #CYPName->ProtIdx
chemfile='/home/hlim/wiZAN/CYP450/PubChem/chemicals/cyp450_ChemIdx.csv'
chemSIDs=[]
sid2idx={}
idx2sid={}
#protnames=['CYP1A2','CYP2C9','CYP2C19','CYP2D6','CYP3A4']
with open(chemfile,"r") as inf:
    header=inf.readline()
    for line in inf:
        line=line.strip().split(",")
        idx=int(line[0])
        sid=str(line[1])
        chemSIDs.append(sid)
        sid2idx[sid]=idx
        idx2sid[idx]=sid

sub=pd.read_csv(extfile,sep=',',header=0,index_col=0,dtype=config.floatX)
sub=np.matrix(sub)
df=pd.DataFrame(sub.T,columns=chemSIDs,dtype=config.floatX)
#print df
store=pd.io.pytables.HDFStore('./ExtendedFC.h5')
store["Fc"]=df
store.close()

sub=pd.read_csv(subfile,sep=',',header=0,index_col=0,dtype=config.floatX)
sub=np.matrix(sub)
df=pd.DataFrame(sub.T,columns=chemSIDs,dtype=config.floatX)
#print df
store=pd.io.pytables.HDFStore('./SubstructFC.h5')
store["Fc"]=df
store.close()

sub=pd.read_csv(pubfile,sep=',',header=0,index_col=0,dtype=config.floatX)
sub=np.matrix(sub)
df=pd.DataFrame(sub.T,columns=chemSIDs,dtype=config.floatX)
#print df
store=pd.io.pytables.HDFStore('./PubChemFC.h5')
store["Fc"]=df
store.close()

