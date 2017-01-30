import numpy as np
import copy

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
    print "nIter: %d, diff: %f\n"%(realIter, diff)
    return r
