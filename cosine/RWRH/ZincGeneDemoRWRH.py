import numpy as np
#import genZincGeneMat2 as genZGM
#import utils.iterRWR as rwr
from scipy import sparse
from scipy.sparse.csr import csr_matrix
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score
from sklearn import metrics
import copy
import traceback
import time
import sys
import os
"""
compute the co-relation between diseases and genes
@author: Yuechang Liu
"""

def  zgDemoRWRH(B, zincCount, geneCount, opt, Test): 
    start_time = time.time()
    c = opt['c']
    maxIter = opt['maxIter']
    tolerance = opt['tolerance']
    lamdaa = opt['lambdaa']
    result = np.zeros(shape=[zincCount,geneCount])
    result = np.matrix(result)
    (nR, nC) = B.shape
    if nR!=nC:
        print "ERROR: not a square matrix!"
        return (result, timeUsed)
    else:
        
        try:
		MZ = B[0:zincCount, 0:zincCount]
		MG = B[zincCount:zincCount+geneCount, zincCount:zincCount+geneCount]
		MZG = B[0:zincCount, zincCount:zincCount+geneCount]
		MGZ = B[zincCount:zincCount+geneCount, 0:zincCount]
		Bw1 = MZG.sum(axis=1)
		Bx = MZ.sum(axis=1)
		M = (1.0*MZ) / Bx
		Bw1[Bw1>0] = 1
		Bw1 = Bw1 * (1.0-lamdaa)
		MZ = np.multiply(M, Bw1)
		By = MG.sum(axis=1)
		MG = (1.0 * MG) / By
		Bw0 = MGZ.sum(axis=1)
		Bw0[Bw0>0] = 1
		Bw0 = Bw0 * (1.0-lamdaa)
		MG = MG * Bw0
		Bw = MZG.sum(axis=1)
		Bz = MGZ.sum(axis=0)
		MZG = lamdaa*MZG / Bw
		MGZ = lamdaa*MGZ / Bz
		MZ[0:zincCount, 0:zincCount] = M[0:zincCount, 0:zincCount]
		B[0:zincCount, 0:zincCount] = MZ
		B[zincCount:zincCount+geneCount, zincCount:zincCount+geneCount] = MG
		B[0:zincCount, zincCount:zincCount+geneCount] = MZG
		B[zincCount:zincCount+geneCount, 0:zincCount] = MGZ
		B[np.isnan(B)] = 0
        except:
                print "error happened during matrix normalization"
                traceback.print_exc()

        for i in range(0,zincCount): #for each chemical (row)
            try:
    #------------------initialize restart vector-----------------
                restart_vector = np.zeros((nR, 1))
                restart_vector[i] = 1
                
    #------------------compute RWR iteratively----------------------- 
                G = csr_matrix(B)
                r = IterRWR(G, c, restart_vector, maxIter, tolerance)
        
    #------------------compute the similarity between tested disease and genes----------------------- 
                gs = r[zincCount : zincCount + geneCount]
#                orderedGeneList = np.argsort(gs[:,0])[::-1]
#                orderedGeneList = orderedGeneList + zincCount
                result[i,0:] = gs.transpose()
            except:
                print "error happens for  zinc %d, try next..."%(i)
                traceback.print_exc()
	timeUsed = time.time() - start_time
        y_score=np.asarray(result).reshape((1,zincCount*geneCount)).flatten()
        y_true=np.asarray(Test).reshape((1,zincCount*geneCount)).flatten()
        print y_score.shape
        print y_true.shape
        AUC=roc_auc_score(y_true,y_score)
        precision, recall, thresholds = precision_recall_curve(y_true,y_score)
        AUPR=metrics.auc(recall,precision)
        return (AUC,AUPR,timeUsed)

def load_data_from_file_csv(dataset, folder):
    with open('/scratch/hansaim.lim/wiZAN/ZINC_data/chem_chem/chem_chem_zinc.txt', "r") as inf:  # the drug similarity file
        drug_sim = [line.strip("\n").split()[0:] for line in inf]
    Csim=np.matrix(drug_sim, dtype=np.float64)
    with open('/scratch/hansaim.lim/wiZAN/ZINC_data/prot_prot/prot_prot_zinc.txt', "r") as inf:  # the target similarity file
        target_sim = [line.strip("\n").split()[0:] for line in inf]
    Tsim=np.matrix(target_sim, dtype=np.float64)
    n=12384	#number of chemicals in ZINC data
    m=3500	#number of proteins in ZINC data
    with open(os.path.join(folder, "train_"+ dataset+".csv"), "r") as inf:
        row=[]
        col=[]
        val=[]
        for line in inf:
            z=line.strip().split(",")
            row.append(int(z[0])-1)
            col.append(int(z[1])-1)
            val.append(int(1))
        train=sparse.coo_matrix((val,(row,col)), shape=(n,m)).todense()
        
    with open(os.path.join(folder, "test_"+ dataset+".csv"), "r") as inf:
        row=[]
        col=[]
        val=[]
        for line in inf:
            z=line.strip().split(",")
            row.append(int(z[0])-1)
            col.append(int(z[1])-1)
            val.append(int(1))
        test=sparse.coo_matrix((val,(row,col)), shape=(n,m)).todense()
    return train, test, Csim, Tsim

def buildGraph(Train,Csim,Tsim):
    (nR,nC)=Train.shape
    Graph=np.zeros(shape=[nR+nC,nR+nC])
    Graph=np.matrix(Graph)
    Graph[0:nR,0:nR]=Csim
    Graph[nR:nR+nC,nR:nR+nC]=Tsim
    Graph[0:nR,nR:nR+nC]=Train
    Graph[nR:nR+nC,0:nR]=Train.transpose()
    return Graph

def IterRWR(A, c, prefer_vec, maxIter, tolerance):
    (nR, nC) = A.shape
    if nR!=nC:
        print 'ERROR: not a square matrix'
        return
    r = prefer_vec;

    realIter = maxIter;
    for i in range(maxIter):
        old_r = copy.deepcopy(r)
        r = (1-c)*A*r + c*prefer_vec
        diff = np.sum(np.absolute(old_r - r));
        if diff < tolerance:
            realIter = i
            break
#    print "nIter: %d, diff: %f\n"%(realIter, diff)
    return r    

def getArguments(Comm):
    """This method checks commands line arguments"""
    tolerance=1e-5
    numLeft=1
    maxIter=200
    if len(Comm)==5:
        inputdir=str(Comm[1])
        dataset=Comm[2]
        reProb=float(Comm[3]) #restart probability
        lamb=float(Comm[4]) #lambda (weights)        
        opt={'c':reProb, 'maxIter':maxIter, 'tolerance':tolerance, 'lambdaa':lamb }
        Train, Test, Csim, Tsim=load_data_from_file_csv(dataset, inputdir)
        initGraph=buildGraph(Train,Csim,Tsim)
        (dC,tC)=Train.shape
        (AUC,AUPR,time)=zgDemoRWRH(initGraph,dC,tC,opt,Test)
        print "Input directory: %s\nAUC: %s, AUPR: %s, time: %s"%(str(inputdir),str(AUC),str(AUPR),str(time))
    else:
        print "example: python zgDemoRWRH.py inputdir restartProb weight"
        exit(1)
    return

if __name__ == '__main__':
    t1=time.time()
    Comm=sys.argv
    getArguments(Comm)
    t2=time.time()
