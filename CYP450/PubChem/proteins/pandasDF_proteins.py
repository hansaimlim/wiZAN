#!/usr/bin/python

import sys
import numpy as np
import pandas as pd
from theano import config

protIdx={'CYP1A2':1,'CYP2C9':2,'CYP2C19':3,'CYP2D6':4,'CYP3A4':5} #CYPName->ProtIdx
chemfile='/home/hlim/wiZAN/CYP450/PubChem/chemicals/cyp450_ChemIdx.csv'
chemSIDs=[]
protnames=['CYP1A2','CYP2C9','CYP2C19','CYP2D6','CYP3A4']

domsim=pd.read_csv('./human_cyp_features_rnacomm.csv',sep=',',header=None,dtype=config.floatX)
domsim=np.matrix(domsim)
df=pd.DataFrame(domsim,columns=protnames,dtype=config.floatX)
store=pd.io.pytables.HDFStore('./DomainSim.h5')
store["Fp"]=df
store.close()

#extdf=pd.read_csv('./ExtFeat_.csv',sep=',',header=None)
seqsim=pd.read_csv('./cyp_blast_sim.csv',sep=',',header=None,dtype=config.floatX)
seqsim=np.matrix(seqsim)
seqdf=pd.DataFrame(seqsim,columns=protnames,dtype=config.floatX)
store=pd.io.pytables.HDFStore('./SequenceSim.h5')
store["Fp"]=seqdf
store.close()

#pd.io.pytables.read_hdf('./SubstructFeat.h5','ChemFeat/SubstructFeat',where=['index<10'])


#y=pd.DataFrame(np.zeros((len(protnames),len(chemSIDs))),columns=chemSIDs,index=protnames)
#store=pd.io.pytables.HDFStore('./CYP_chemical_interactions.h5')
#store["matrix"]=y
#store.close()

