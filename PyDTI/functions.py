
import os
import numpy as np
from scipy import sparse
import scipy.io as sio
from collections import defaultdict

def get_rcrs(sarr, testPair):
    #smat=score matrix
    #testPair: test pairs passed by zip(row,col) containing test pairs
    #testMat format: drugs in rows; targets in cols (transposed in load_data_from_file function in this class
    #rank is dense rank on each row
    smat=np.asmatrix(sarr).reshape((12384,3500))
    rcrs=[]
    for d, t in testPair:
        boo=smat[d,:]>=smat[d,t]
        rank=np.sum(boo,dtype=np.int32)
        r=[d,t,rank,smat[d,t]]
        rcrs.append(r)
    return rcrs
def TPR_by_cutRank(rcrs, cutRank):
    #True positive rate at cutoff rank
    #rcrs=[row,col,rank,score]
    #cutRank=cutoff Rank (integer)

    tp=0
    for ll in rcrs:
        if ll[2] <= cutRank:
            tp+=1
    tpr= float(tp)/float(len(rcrs))
    return tpr

def load_data_from_file_demo(dataset, folder):
    testdd=sio.loadmat('/scratch/hansaim.lim/data/testdat/testdd.mat')['testdd']
    testtt=sio.loadmat('/scratch/hansaim.lim/data/testdat/testtt.mat')['testtt']
    n=5	#number of chemicals in demo data
    m=4	#number of proteins in demo data
    row=[]
    col=[]
    val=[]
    with open(os.path.join(folder, "train_"+ dataset + ".csv"), "r") as tr:
        for line in tr:
            z=line.strip().split(",")
            row.append(int(z[1])-1)	#0-based index
            col.append(int(z[0])-1)
            val.append(1)
    train=sparse.coo_matrix((val,(row,col)), shape=(m, n)).toarray()
    trainMat=np.array(train, dtype=np.float64).T
    
    row=[]
    col=[]
    val=[]
    with open(os.path.join(folder, "test_"+ dataset + ".csv"), "r") as ts:
        for line in ts:    
            z=line.strip().split(",")
            row.append(int(z[1])-1) #0-based index
            col.append(int(z[0])-1)
            val.append(1)
    test=sparse.coo_matrix((val,(row,col)), shape=(m, n)).toarray()
    testMat=np.array(test, dtype=np.float64).T
    return trainMat, testMat, testdd, testtt

def load_data_from_file_csv(dataset, folder):
    with open('/scratch/hansaim.lim/wiZAN/ZINC_data/chem_chem/chem_chem_zinc.txt', "r") as inf:  # the drug similarity file
        drug_sim = [line.strip("\n").split()[0:] for line in inf]

    with open('/scratch/hansaim.lim/wiZAN/ZINC_data/prot_prot/prot_prot_zinc.txt', "r") as inf:  # the target similarity file
        target_sim = [line.strip("\n").split()[0:] for line in inf]
    n=12384	#number of chemicals in ZINC data
    m=3500	#number of proteins in ZINC data
    with open(os.path.join(folder, "train_"+ dataset+".csv"), "r") as inf:
        row=[]
        col=[]
        val=[]
        for line in inf:
            z=line.strip().split(",")
            row.append(int(z[1])-1)
            col.append(int(z[0])-1)
            val.append(1)
        int_array=sparse.coo_matrix((val,(row,col)), shape=(m,n)).toarray()
        
    with open(os.path.join(folder, "test_"+ dataset+".csv"), "r") as inf:
        row=[]
        col=[]
        val=[]
        for line in inf:
            z=line.strip().split(",")
            row.append(int(z[1])-1)
            col.append(int(z[0])-1)
            val.append(1)
        test_array=sparse.coo_matrix((val,(row,col)), shape=(m,n)).toarray()

    intMat = np.array(int_array, dtype=np.float64).T    # drug-target interaction matrix
    testMat = np.array(test_array, dtype=np.float64).T    # drug-target interaction matrix
    drugMat = np.array(drug_sim, dtype=np.float64)      # drug similarity matrix
    targetMat = np.array(target_sim, dtype=np.float64)  # target similarity matrix
    return intMat, testMat, drugMat, targetMat

def load_data_from_file(dataset, folder):
    with open(os.path.join(folder, dataset+"_admat_dgc.txt"), "r") as inf:
        inf.next()
        int_array = [line.strip("\n").split()[1:] for line in inf]

    drug_sim=sio.loadmat('/scratch/hansaim.lim/data/testdat/testdd.mat')['testdd']
    target_sim=sio.loadmat('/scratch/hansaim.lim/data/testdat/testtt.mat')['testtt']
#    with open(os.path.join(folder, dataset+"_simmat_dc.txt"), "r") as inf:  # the drug similarity file
#        inf.next()
#        drug_sim = [line.strip("\n").split()[1:] for line in inf]
#
#    with open(os.path.join(folder, dataset+"_simmat_dg.txt"), "r") as inf:  # the target similarity file
#        inf.next()
#        target_sim = [line.strip("\n").split()[1:] for line in inf]

    intMat = np.array(int_array, dtype=np.float64).T    # drug-target interaction matrix
    drugMat = np.array(drug_sim, dtype=np.float64)      # drug similarity matrix
    targetMat = np.array(target_sim, dtype=np.float64)  # target similarity matrix
    return intMat, drugMat, targetMat

def get_drug_target_names_demo(*args):
    drugname=[]
    targetname=[]
    for dline in open("/scratch/hansaim.lim/data/testdat/demo_dlist.tsv","r").xreadlines():
        d=dline.strip().split("\t")
        dind=int(d[0])
        dname=str(d[1])
        drugname.append(dname)
    for tline in open("/scratch/hansaim.lim/data/testdat/demo_tlist.tsv","r").xreadlines():
        t=tline.strip().split("\t")
        tind=int(t[0])
        tname=str(t[1])
        targetname.append(tname)
    return drugname, targetname

def get_drug_target_names_zinc(*args):
    drugname=[]
    targetname=[]
    for dline in open("/scratch/hansaim.lim/wiZAN/ZINC_data/ZINC_chemical.tsv","r").xreadlines():
        d=dline.strip().split("\t")
        dind=int(d[0])
        dname=str(d[1])
        drugname.append(dname)
    for tline in open("/scratch/hansaim.lim/wiZAN/ZINC_data/ZINC_protein_index.tsv","r").xreadlines():
        t=tline.strip().split("\t")
        tind=int(t[0])
        tname=str(t[1])
        targetname.append(tname)
    return drugname, targetname
    
def get_drugs_targets_names(dataset, folder):
    with open(os.path.join(folder, dataset+"_admat_dgc.txt"), "r") as inf:
        drugs = inf.next().strip("\n").split()
        targets = [line.strip("\n").split()[0] for line in inf]
    return drugs, targets


def cross_validation(intMat, seeds, cv=0, num=10):
    cv_data = defaultdict(list)
    for seed in seeds:
        num_drugs, num_targets = intMat.shape
        prng = np.random.RandomState(seed)
        if cv == 0:
            index = prng.permutation(num_drugs)
        if cv == 1:
            index = prng.permutation(intMat.size)
        step = index.size/num
        for i in xrange(num):
            if i < num-1:
                ii = index[i*step:(i+1)*step]
            else:
                ii = index[i*step:]
            if cv == 0:
                test_data = np.array([[k, j] for k in ii for j in xrange(num_targets)], dtype=np.int32)
            elif cv == 1:
                test_data = np.array([[k/num_targets, k % num_targets] for k in ii], dtype=np.int32)
            x, y = test_data[:, 0], test_data[:, 1]
            test_label = intMat[x, y]
            W = np.ones(intMat.shape)
            W[x, y] = 0
            cv_data[seed].append((W, test_data, test_label))
    return cv_data


def train(model, cv_data, intMat, drugMat, targetMat):
    aupr, auc = [], []
    for seed in cv_data.keys():
        for W, test_data, test_label in cv_data[seed]:
            model.fix_model(W, intMat, drugMat, targetMat, seed)
            aupr_val, auc_val = model.evaluation(test_data, test_label)
            aupr.append(aupr_val)
            auc.append(auc_val)
    return np.array(aupr, dtype=np.float64), np.array(auc, dtype=np.float64)


def svd_init(M, num_factors):
    from scipy.linalg import svd
    U, s, V = svd(M, full_matrices=False)
    ii = np.argsort(s)[::-1][:num_factors]
    s1 = np.sqrt(np.diag(s[ii]))
    U0, V0 = U[:, ii].dot(s1), s1.dot(V[ii, :])
    return U0, V0.T


def mean_confidence_interval(data, confidence=0.95):
    import scipy as sp
    import scipy.stats
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, h


def write_metric_vector_to_file(auc_vec, file_name):
    np.savetxt(file_name, auc_vec, fmt='%.6f')


def load_metric_vector(file_name):
    return np.loadtxt(file_name, dtype=np.float64)
