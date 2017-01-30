import numpy as np
from scipy import sparse
from scipy.sparse.csr import csr_matrix
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score
from sklearn import metrics
import copy
import traceback
import sys
import os

from timeit import default_timer as timer

from sklearn.cross_validation import train_test_split
from sklearn.cross_validation import KFold
"""
compute the co-relation between diseases and genes
@author: Yuechang Liu
"""

def  RWRH(B, zincCount, geneCount, opt, Test, TestChemIdx): 
    c = opt['c']
    maxIter = opt['maxIter']
    tolerance = opt['tolerance']
    lamdaa = opt['lambdaa']
    result = np.zeros(shape=Test.shape)
    result = np.matrix(result)
    print Test.shape
    print result.shape
    (nR, nC) = B.shape
    if nR!=nC:
        print "ERROR: not a square matrix!"
        return result
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

        for ii in range(0,len(TestChemIdx)): #for each test chemical
            chemIdx=TestChemIdx[ii]
            try:
    #------------------initialize restart vector-----------------
                restart_vector = np.zeros((nR, 1))
                restart_vector[chemIdx] = 1
                
    #------------------compute RWR iteratively----------------------- 
                G = csr_matrix(B)
                r = IterRWR(G, c, restart_vector, maxIter, tolerance)
        
    #------------------compute the similarity between tested disease and genes----------------------- 
                gs = r[zincCount : zincCount + geneCount]
#                orderedGeneList = np.argsort(gs[:,0])[::-1]
#                orderedGeneList = orderedGeneList + zincCount
                result[ii,0:] = gs.transpose()
            except:
                print "error happens for  zinc %d, try next..."%(ii)
                traceback.print_exc()
        y_score=np.asarray(result).reshape(len(TestChemIdx)*geneCount).flatten()
        y_true=np.asarray(Test).reshape(len(TestChemIdx)*geneCount).flatten()
        AUC=roc_auc_score(y_true,y_score)
        precision, recall, thresholds = precision_recall_curve(y_true,y_score)
        AUPR=metrics.auc(recall,precision)
        return (AUC,AUPR)

def load_matrix(TargetTargetSimFile,DrugDrugSimFile,TargetDrugNetworkFile):
    with open(DrugDrugSimFile, "r") as inf:  # the drug similarity file
        dHead=inf.readline()
        numDrug=0
        drug_sim=[]
        for line in inf:
            line=line.rstrip('\r\n')
            if not line:
                continue
            row=line.strip().split("\t")
            drug_sim.append([float(r) for r in row[1:]])
            numDrug+=1
    Csim=np.array(drug_sim, dtype=np.float64).reshape(numDrug,numDrug)
    Csim=np.matrix(Csim)
    with open(TargetTargetSimFile, "r") as inf:  # the target similarity file
        tHead=inf.readline()
        numTarg=0
        target_sim=[]
        for line in inf:
            line=line.rstrip('\r\n')
            if not line:
                continue
            row=line.strip().split("\t")
            target_sim.append([float(r) for r in row[1:]])
            numTarg+=1
    Tsim=np.array(target_sim, dtype=np.float64).reshape(numTarg,numTarg)
    Tsim=np.matrix(Tsim)
    with open(TargetDrugNetworkFile, "r") as inf:  # the target-drug adjacency matrix
        tdHead=inf.readline()
        target_drug=[]
        for line in inf:
            line=line.rstrip('\r\n')
            if not line:
                continue
            row=line.strip().split("\t")
            target_drug.append([float(r) for r in row[1:]])
    Train=np.array(target_drug, dtype=np.float64).reshape(numTarg,numDrug)
    Train=np.matrix(Train).transpose()
    return Train, Csim, Tsim

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
    if len(Comm)==5:
        start = timer()
        ttsim=str(Comm[1])
        ddsim=str(Comm[2])
        tdnet=str(Comm[3])
        Nfold=int(Comm[4])
        Whole, Csim, Tsim=load_matrix(ttsim,ddsim,tdnet)
        (numChem,numTarg)=Whole.shape
        opt={'c':0.7, 'maxIter':100, 'tolerance':1e-9, 'lambdaa':0.5 }

        #find non-zero pairs
        row,col=np.where(Whole>0)
        pairs=[]
        for i in range(0,row.shape[1]):
            pairs.append([row[0,i],col[0,i]])
        
        #Set N-fold iterator
        kf = KFold(len(pairs), n_folds=Nfold, shuffle=True)
        AUC=[]
        AUPR=[]
        for train_pair_index, test_pair_index in kf:
            #run N-fold validation
            Train=copy.deepcopy(Whole)
            Test=None
            for idx in test_pair_index: #masking each test pair
                [testrow,testcol]=pairs[idx]
                Train[testrow,testcol]=0 #masked
            Test=Whole-Train
            testrows, testcols = np.where(Test>0)
            TestChemIdx=list(set(testrows))
            TestChemIdx=list(TestChemIdx[0])
            TestCompressed = Test[TestChemIdx]

            initGraph=buildGraph(Train,Csim,Tsim)
            (auc,aupr)=RWRH(initGraph,numChem,numTarg,opt,TestCompressed,TestChemIdx)
            AUC.append(auc)
            AUPR.append(aupr)

        end=timer()
        timeUsed=end - start
        AUCmean=np.mean(AUC)
        AUCstd=np.std(AUC)
        AUPRmean=np.mean(AUPR) 
        AUPRstd=np.std(AUPR)       
        print "Target-Target:%s,  Drug-Drug:%s,  Target-Drug:%s,  K-fold:%s"%(str(ttsim),str(ddsim),str(tdnet),str(Nfold))
        print "mean AUC: %s (std. %s), mean AUPR: %s (std. %s), time: %s"%(str(AUCmean),str(AUCstd),str(AUPRmean),str(AUPRstd),str(timeUsed))
        print "AUC values:"
        print AUC
    else:
        print "example: python zgDemoRWRH.py target-target.txt drug-drug.txt target-drug.txt 5 (N-fold)"
        exit(1)
    return

if __name__ == '__main__':
    Comm=sys.argv
    getArguments(Comm)
