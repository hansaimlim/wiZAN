#!/usr/bin/python
# encoding: utf-8
#modified by Hansaim Lim hansaimlim@gmail.com on 5/27/2015
#calculate TPR by cutoff row rank
#print out TPR by cutoff row rank, instead of predicted scores for all pairs

# Written by Alexios Koutsoukas ak726@cam.ac.uk
#Supervisors: Profesor RC Glen & Dr. A. Bender
#Funding: Unilever
#all rights reserved 2012
#Instructions: to run this program: python PRW.py TrainFile TestFile OutputFile h
#where h the smoothing parameter [0.1,0.99]
#Training and testfile needs to be in tab delimited fromat: 
#trainfile first column must contain the classes ids
#test file should be in same format with trainfile
#This software is provided with publication:
# Koutsoukas et al. In silico target predictions: defining a benchmarking dataset and comparison
#of performance of the multiclass Naïve Bayes and Parzen-Rosenblatt Window. JCIM 2013

#libraries to be imported
import random
import time
import os
import sys
import math
import numpy as np
import sets
#library for TPR calculation
from collections import defaultdict

def get_feature():
    #chemical features
    Feature={}
    feature='/home/hlim/wizan/wiZAN/ZINC_data/zinc_chemIndex_feature.tsv'
    try:
        for line in open(feature,"r").xreadlines():
            f=line.strip().split("\t")
            chem=str(f[0])
            feat=f[1:]
            Feature[chem]=feat
    except:
        print "Please provide a vaild chemical feature file"
        exit(1)
    return Feature
def get_train(trainfile,Feature):
    #Molecules={} training instances
    Train={}
    #Number of instances
    count=0
    print "Reading training instances ",trainfile
    try:
        for line in open(trainfile,"r").xreadlines():
            z=line.strip().split(",")
            count+=1
#            if count%5000==0:
#                print count
            chem=str(z[0])
            prot=str(z[1])
            feat=Feature[chem]
            Train.setdefault(prot,[]).append(feat)
    except:
        print "Please provide a valid file"
        exit(1)
#    print "Total number of training instances: ",count
#    print "Number of classes: ",len(Train.keys())
    return Train


def DifferenceVectors(alist,blist):
    """Calculates difference of vectors
    and returns d"""
    #alist a molecule
    #blist b molecule
    #d difference of two sparce vectors in power of 2
    A=set(alist)
    B=set(blist)
    D=int(len(A.difference(B))+len(B.difference(A)))
    return D

def GaussianKernel(d,l):
    """Gaussian kernel:
    l=smoothing parameter
    """
    K=float(math.exp(-(float(d)/float(2*math.pow(l,2)))))
    return K

def get_PRW_rcrs(Feature,Train,testPair,l):
    #smat=score matrix
    #testPair: test pairs passed by zip(row,col) containing test pairs
    #testMat format: drugs in rows; targets in cols (transposed in load_data_from_file function in this class
    #rank is dense rank on each row
    rcrs=[]
    for row, col in testPair:
        row=str(row)
        col=str(col)
        predscore=0.0
        rowrank=1
        Pred={} #dictionary holding the predictions for each class
        feat=Feature[row]
        sumTotal=0
        for key in Train:
            Dist=[] #will hold kernelized distances to training class
            for j in Train[key]:
                d=DifferenceVectors(feat,j)
                K=GaussianKernel(d,l)
                Dist.append(float(K))
            sumTotal+=float(sum(Dist))
            Pred[key]=float(sum(Dist))
        ratioPred={} # dictionary holding the normalized propabilities for each class
        for key in Pred:
            ratioPred[key]=float(Pred[key])/float(sumTotal)
        for i in sorted(ratioPred.iteritems(),key=lambda x:x[1],reverse=True):  #find predicted score for test pair
            if str(i[0])==col:
	        predscore=float(i[1])
	        break
        for i in sorted(ratioPred.iteritems(),key=lambda x:x[1],reverse=True):  #count dense rank
            if float(i[1])<predscore:
	        #score smaller than predict score. does not increase the rank
	        continue
            else:      #rank decrease (equal or higher #) if a score found greater or equal to the predicted score (Dense Rank)
                rowrank+=1
        r=[row, col, rowrank, predscore]
        rcrs.append(r)
    return rcrs
def TPRbyRowRank(rcrs, cutRank):
    #True positive rate at cutoff rank
    #rcrs=[row,col,rank,score]
    #cutRank=cutoff Rank (integer)

    tp=0
    for ll in rcrs:
        if ll[2] <= cutRank:
            tp+=1
    tpr= float(tp)/float(len(rcrs))
    return tpr

def getArguments(Comm):
    """This method checks commands line arguments"""
    input_dir=''
    dataset=''
    l=0.0
    if len(Comm)==4:
        input_dir=str(Comm[1])
        dataset=str(Comm[2])
        l=float(Comm[3])
    else:
        print "To run this program you have to provide a training file, a test file, an output file and a smoothing parameter"
        print "the chemicals must be integer indexed. This script read chemical feature files based on chemical index in integer"
        print "example: python PRW_10cv.py /path/to/N1/L1to5/ N1L1to5 0.9"
        exit(1)
    feature=get_feature()
    RCRS=[]
    TPR35=[]
    for i in range(1,11):
        #10cv for 1 to 10
        trainfile=input_dir+"train"+str(i)+".csv"
        Train=get_train(trainfile,feature)

        testfile=input_dir+"test"+str(i)+".csv"
         
        testrow=[]
        testcol=[]
        print "Reading test instances : ",testfile
        for line in open(testfile,"r").xreadlines():
            z=line.strip().split(",")
            chem=str(z[0])
            prot=str(z[1])
            testrow.append(chem)
            testcol.append(prot)
        rcrs=get_PRW_rcrs(feature,Train,zip(testrow,testcol),l)
        tpr35=TPRbyRowRank(rcrs,35)
        TPR35.append(tpr35)
        RCRS=RCRS+rcrs
    avgtpr35=np.average(TPR35)
    semtpr35=(np.std(TPR35)/math.sqrt(len(TPR35)))
    print "Dataset=%s, Avg.TPR35=%s, S.E.M.TPR35=%s"%(dataset,str(avgtpr35),str(semtpr35))
    print TPR35
    print "Rank\tTPR"
    for i in range(1,351):
        tpr=TPRbyRowRank(RCRS,i)
        print "%s\t%s"%(str(i),str(tpr))
    return
if __name__ == '__main__':
    t1=time.time()
    Comm=sys.argv
    getArguments(Comm)
    t2=time.time()

