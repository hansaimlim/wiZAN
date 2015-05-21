#!/usr/bin/python
# encoding: utf-8

#Input modified by Hansaim Lim on 5/21/2015
#No more test file is needed as an input. Chemical feature file is used instead.
#This script runs PRW test on all unique chemicals in (train + test) against train data
#A set of unique chemicals obtained from the chemical feature file.

#Output modified by Hansaim Lim on 5/20/2015
#csv output for each prediction. (label, target, probability) for each line.

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


def PRW(train,feature,Output,l):
    #chemical features
    Feature={}
    #Molecules={} training instances
    Train={}
    #Number of instances
    count=0
    print "Reading training instances ",train
    try:
	for line in open(feature,"r").xreadlines():
            f=line.strip().split("\t")
            chem=f[0]
            feat=f[1:]
            Feature[chem]=feat
    except:
        print "Please provide a vaild chemical feature file"
        exit(1)
    try:
        for line in open(train,"r").xreadlines():
            #train file in csv format (chem_index, prot_index, 1)
            z=line.strip().split(", ")
            count+=1
            if count%5000==0:
                print count
            chem=z[0]	#chemical index
            label=z[1]	#protein index
            feat=Feature[chem]
            Train.setdefault(label,[]).append(feat)
    except:
        print "Please provide a valid file"
        exit(1)
    print "Done reading training instances"
    print "Total number of training instances: ",count
    print "Number of classes: ",len(Train.keys())
    print "===================================================="
    print "Perform Classification on test instances"
    S=open(Output,"w")
    cow=0
    print "Writing Results to : ",Output
    for label in Feature:
        #all unique chemicals are tested
        cow+=1
        if cow%100==0:
            print cow
        Pred={} #dictionary holding the predictions for each class
        feat=Feature[label]
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
        string=""
        for i in sorted(ratioPred.iteritems(),key=lambda x:x[1],reverse=True):
             S.write(str(label)+", "+str(i[0])+", "+str(i[1])+"\n")
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
    """Method updated by Hansaim Lim on 5/20/2015"""
    """Takes third argument for Chemical Feature file"""
    if len(Comm)==5:
        Train=Comm[1]
	Feature=Comm[2]
        Output=Comm[3]
        l=float(Comm[4])
        PRW(Train,Feature,Output,l)
    else:
        print "To run this program you have to provide a training file, a test file, a feature file, an output file and a smoothing parameter l"
        print "example: python prw.py train test feature output1 0.9"
        exit(1)
    return
    

if __name__ == '__main__':
    t1=time.time()
    Comm=sys.argv
    getArguments(Comm)
    t2=time.time()
