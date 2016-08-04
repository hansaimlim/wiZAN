import numpy as np
from scipy.sparse import coo_matrix, 
from sklearn.metrics import precision_recall_curve, roc_curve
from sklearn.metrics import auc


class WIZAN:

    def __init__(self, reg=0.1, w=0.1, p=0.01, rank=300, max_iter=400, cimp=0.75, pimp=0.1):
        self.reg = float(reg)  
        self.w = float(w)
        self.p = float(p)
        self.rank = int(rank)
        self.max_iter = int(max_iter)
        self.cimp = float(cimp)
        self.pimp = float(pimp)


    def fix_model(self, W, intMat, drugMat, targetMat, seed=None):
        self.num_drugs, self.num_targets = intMat.shape

    def updateUV(self, intMat, Lu, Lv):

        return U, V

    def get_UVT(self, intMat, U, V):
    
        return UVT

    def updateU(self, intMat, UVT, w, p, Lu_plus, Lu_minus, U0, V, imp1, imp2):

        U1(np.isnan(U1)) = 0
        return U1
    def predict_scores(self, test_data, N):
        dinx = np.array(list(self.train_drugs))
        DS = self.dsMat[:, dinx]
       # print DS drug-drug sim with 0 diagonal entries
        tinx = np.array(list(self.train_targets))
        TS = self.tsMat[:, tinx]
       # print TS target-target sim with 0 diagonal entries
        scores = []
        for d, t in test_data:
            if d in self.train_drugs: 
                if t in self.train_targets:
                    val = np.sum(self.U[d, :]*self.V[t, :])
                else:
                    jj = np.argsort(TS[t, :])[::-1][:N]
                    val = np.sum(self.U[d, :]*np.dot(TS[t, jj], self.V[tinx[jj], :]))/np.sum(TS[t, jj])
            else:
                if t in self.train_targets:
                    ii = np.argsort(DS[d, :])[::-1][:N]
                    val = np.sum(np.dot(DS[d, ii], self.U[dinx[ii], :])*self.V[t, :])/np.sum(DS[d, ii])
                else:
                    ii = np.argsort(DS[d, :])[::-1][:N]
                    jj = np.argsort(TS[t, :])[::-1][:N]
                    v1 = DS[d, ii].dot(self.U[dinx[ii], :])/np.sum(DS[d, ii])
                    v2 = TS[t, jj].dot(self.V[tinx[jj], :])/np.sum(TS[t, jj])
                    val = np.sum(v1*v2)
            if np.isnan(val):
                scores.append(0)
            else:
                scores.append(np.exp(val)/(1+np.exp(val)))
       # print smat #whole prediction matrix
        return np.array(scores) #from original code

    def __str__(self):
        return "Model: wizan, reg:%s, wt:%s, imp:%s, rank:%s, p_chem:%s, p_prot:%s, max_iter:%s" % (self.reg, self.w, self.p, self.rank, self.cimp, self.pimp, self.max_iter)
