import pickle
import scipy.sparse
import numpy as np
import itertools
from sklearn.decomposition import NMF
import time as time


def matrix_factorization(R, P, Q, K, steps=5000, alpha=0.0002, beta=0.02):
    Q = Q.T
    for step in xrange(steps):
        for i in xrange(len(R)):
            for j in xrange(len(R[i])):
                if R[i][j] > 0:
                    eij = R[i][j] - np.dot(P[i,:],Q[:,j])
                    for k in xrange(K):
                        P[i][k] = P[i][k] + alpha * (2 * eij * Q[k][j] - beta * P[i][k])
                        Q[k][j] = Q[k][j] + alpha * (2 * eij * P[i][k] - beta * Q[k][j])
        eR = np.dot(P,Q)
        e = 0
        for i in xrange(len(R)):
            for j in xrange(len(R[i])):
                if R[i][j] > 0:
                    e = e + pow(R[i][j] - np.dot(P[i,:],Q[:,j]), 2)
                    for k in xrange(K):
                        e = e + (beta/2) * (pow(P[i][k],2) + pow(Q[k][j],2))
        if e < 0.001:
            break
    return P, Q.T


filepath='../Zeisel_pipeline/'

with open(filepath+'Zeisel_TCC_pairwise_JS_distance_subsample100_full.dat','rb') as infile:
    D = pickle.load(infile)

# Sanity check
assert np.all(np.isclose(D,D.T))
assert np.all(np.isclose(np.diag(D),np.zeros(np.diag(D).shape)))

# Zeisel's 9 main clusters
labels9 = np.loadtxt(filepath+'./Zeisels_data/Zeisels_labels9.txt',dtype=str).astype(int)-1

# Zeisel's neurons (labeled as 0) and non-neurons (labeled as 1) obtained from labels9
labels2 = np.copy(labels9)
for i in [0,1,2]: labels2[labels2 == i] = 0
labels2[labels2 != 0] = 1

# Zeisel's 47 total clusters
labels47 = np.loadtxt(filepath+'./Zeisels_data/Zeisels_labels47.txt',dtype=str)    

a = time.time()
model = NMF(n_components=5, init='random', random_state=0)

W = model.fit_transform(D)
H = model.components_
print time.time() - a

D_2 = np.dot(W, H)
diff = D - D_2
print(diff)
print(D_2.shape)

# N = len(D)
# M = len(D[0])
# K = 2

# P = np.random.rand(N,K)
# Q = np.random.rand(M,K)

# nP, nQ = matrix_factorization(D, P, Q, K)
# nR = np.dot(nP, nQ.T)
# print(nR)

