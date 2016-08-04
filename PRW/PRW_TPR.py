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
import sets
#library for TPR calculation
from collections import defaultdict


def PRW(train,test,Output,l):
    #chemical features
    Feature={}
    feature='/home/hansaim/wiZAN/ZINC_data/zinc_chemIndex_feature.tsv'
    try:
        for line in open(feature,"r").xreadlines():
            f=line.strip().split("\t")
            chem=f[0]
            feat=f[1:]
            Feature[chem]=feat
    except:
        print "Please provide a vaild chemical feature file"
        exit(1)
    #Molecules={} training instances
    Train={}
    #Number of instances
    count=0
    print "Reading training instances ",train
    try:
        for line in open(train,"r").xreadlines():
            z=line.strip().split(", ")
            count+=1
            if count%5000==0:
                print count
            chem=z[0]
            prot=z[1]
            feat=Feature[chem]
            Train.setdefault(prot,[]).append(feat)
    except:
        print "Please provide a valid file"
        exit(1)
    print "Done reading training instances"
    print "Total number of training instances: ",count
    print "Number of classes: ",len(Train.keys())
    print "===================================================="
    print "Perform Classification on test instances"
    print "Reading test instances : ",test
    S=open(Output,"w")
    cp=0	#condition positive
    print "Writing Results to : ",Output
    testRank={}	#dictionary holding the test pair and its row rank (dense rank)
    for line in open(test,"r").xreadlines():
        cp+=1
        if cp%100==0:
            print cp
        z=line.strip().split(", ")
        Pred={} #dictionary holding the predictions for each class
        chem=z[0]
        prot=str(z[1])
        predscore=float(0.0)
        rowrank=1 #rank starts from 1 and increases
        feat=Feature[chem]
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
        for i in sorted(ratioPred.iteritems(),key=lambda x:x[1],reverse=True):	#find predicted score for test pair
             if str(i[0])==prot:
                  predscore=float(i[1])
                  break
        for i in sorted(ratioPred.iteritems(),key=lambda x:x[1],reverse=True):	#count dense rank
             if float(i[1])<predscore:
                  #score smaller than predict score. does not increase the rank
                  continue
             else:	#rank decrease (equal or higher #) if a score found greater or equal to the predicted score (Dense Rank)
                  rowrank+=1
        pair=str(chem)+", "+str(prot)
        testRank[pair]=rowrank	#assign rank for the test pair
    rankTP=defaultdict(int)	#dictionary for True Positive counts by cutoff rowrank
    for chem_prot in testRank:
        if int(testRank[chem_prot])>100:	#no need to count if ranked over 100
             continue
        for cutoff in range(1, 101):	#count tp for up to rank 100
             if int(testRank[chem_prot])<=int(cutoff):
                  rankTP[cutoff]+=1
    for cut in rankTP:	#print out TPR by rowRank
        tpr=float(rankTP[cut])/float(cp)	#true positive / condition positive
        S.write(str(cut)+"\t"+str(tpr)+"\n")
    return 
    

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

def getArguments(Comm):
    """This method checks commands line arguments"""
    if len(Comm)==5:
        Train=Comm[1]
        Test=Comm[2]
#        Feature=Comm[3]
        Output=Comm[3]
        l=float(Comm[4])
        PRW(Train,Test,Output,l)
    else:
        print "To run this program you have to provide a training file, a test file, an output file and a smoothing parameter"
	print "the chemicals must be integer indexed. This script read chemical feature files based on chemical index in integer"
        print "example: python prw.py train test output1 0.9"
        exit(1)
    return
    

if __name__ == '__main__':
    t1=time.time()
    Comm=sys.argv
    getArguments(Comm)
    t2=time.time()
    

