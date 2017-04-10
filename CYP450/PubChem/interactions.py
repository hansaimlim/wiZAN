#!/usr/bin/python

import sys
import numpy as np
import pandas as pd

def readCYP(filename,chemIdx):
    #the argument file lists chemical assay data for a given CYP450 protein target
    chemIdxs=[] #list of chemical index interacting with the given protein (one protein per file)
    with open(filename, "r") as inf:
        header=inf.readline()
        #print header
        for line in inf:
            line=line.strip().split(",")
            sid=str(line[1])
            idx=chemIdx[sid]
            chemIdxs.append(sid)
    return chemIdxs

files=['./cyp1a2.csv', './cyp2c19.csv',  './cyp2c9.csv',  './cyp2d6.csv',  './cyp3a4.csv']

chemIdx={} #PubchemSID->ChemIdx
protIdx={'CYP1A2':1,'CYP2C9':2,'CYP2C19':3,'CYP2D6':4,'CYP3A4':5} #CYPName->ProtIdx
chemfile='/home/hlim/wiZAN/CYP450/PubChem/chemicals/cyp450_ChemIdx.csv'
chemSIDs=[]
protnames=['CYP1A2','CYP2C9','CYP2C19','CYP2D6','CYP3A4']
with open(chemfile, "r") as inf:
    header=inf.readline()
    for line in inf:
        line=line.strip().split(',')
        idx=int(line[0])
        sid=str(line[1])
        chemIdx[sid]=idx
        chemSIDs.append(sid)
#print len(protnames) #5
#print len(chemSIDs) #17143

y=pd.DataFrame(np.zeros((len(protnames),len(chemSIDs))),columns=chemSIDs,index=protnames)
CYP1A2=readCYP('./cyp1a2.csv',chemIdx)
for c in CYP1A2:
    y[c]['CYP1A2']=1
CYP2C9=readCYP('./cyp2c9.csv',chemIdx)
for c in CYP2C9:
    y[c]['CYP2C9']=1
CYP2C19=readCYP('./cyp2c19.csv',chemIdx)
for c in CYP2C19:
    y[c]['CYP2C19']=1
CYP2D6=readCYP('./cyp2d6.csv',chemIdx)
for c in CYP2D6:
    y[c]['CYP2D6']=1
CYP3A4=readCYP('./cyp3a4.csv',chemIdx)
for c in CYP3A4:
    y[c]['CYP3A4']=1
store=pd.io.pytables.HDFStore('./CYP_chemical_interactions.h5')
store["matrix"]=y
store.close()

