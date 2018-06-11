import pickle
import scipy.sparse
import numpy as np
import itertools
from sklearn.decomposition import NMF
import time as time
from sklearn import manifold
from sklearn import cluster
import networkx as nx


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

def tSNE_pairwise(D):
    tsne = manifold.TSNE(n_components=2, random_state=0, metric='precomputed', n_iter=2000, verbose=1);
    X_tsne = tsne.fit_transform(D);
    return X_tsne
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


# Spectral clustering with TCCs achieved an error rate of 2.59567387687% (2 labels).
# Spectral clustering with TCCs achieved an error rate of 38.7021630616% (9 labels).


# model = NMF(n_components=5, init='random', random_state=0)
# Spectral clustering with TCCs achieved an error rate of 2.56239600666% (2_2 labels).
# Spectral clustering with TCCs achieved an error rate of 50.1164725458% (9_2 labels).

# model = NMF(n_components=10, init='random', random_state=0)
# Spectral clustering with TCCs achieved an error rate of 2.62895174709% (2_2 labels).
# Spectral clustering with TCCs achieved an error rate of 47.7204658902% (9_2 labels).

model = NMF(n_components=13, init='random', random_state=0)

# model = NMF(n_components=14, init='random', random_state=0)
# Spectral clustering with TCCs achieved an error rate of 2.59567387687% (2_2 labels).
# Spectral clustering with TCCs achieved an error rate of 36.8053244592% (9_2 labels).

# model = NMF(n_components=15, init='random', random_state=0)
# Spectral clustering with TCCs achieved an error rate of 2.59567387687% (2_2 labels).
# Spectral clustering with TCCs achieved an error rate of 36.3727121464% (9_2 labels).

# model = NMF(n_components=16, init='random', random_state=0)
# Spectral clustering with TCCs achieved an error rate of 2.59567387687% (2_2 labels).
# Spectral clustering with TCCs achieved an error rate of 43.0615640599% (9_2 labels).

# model = NMF(n_components=17, init='random', random_state=0)
# Spectral clustering with TCCs achieved an error rate of 2.59567387687% (2_2 labels).
# Spectral clustering with TCCs achieved an error rate of 39.4675540765% (9_2 labels).

# model = NMF(n_components=18, init='random', random_state=0)
# Spectral clustering with TCCs achieved an error rate of 2.59567387687% (2_2 labels).
# Spectral clustering with TCCs achieved an error rate of 39.3011647255% (9_2 labels).
# model = NMF(n_components=19, init='random', random_state=0)

# model = NMF(n_components=20, init='random', random_state=0)
# Spectral clustering with TCCs achieved an error rate of 2.59567387687% (2_2 labels).
# Spectral clustering with TCCs achieved an error rate of 39.1680532446% (9_2 labels).

# model = NMF(n_components=25, init='random', random_state=0)
# Spectral clustering with TCCs achieved an error rate of 2.59567387687% (2_2 labels).
# Spectral clustering with TCCs achieved an error rate of 38.8019966722% (9_2 labels).

# model = NMF(n_components=50, init='random', random_state=0)
# Spectral clustering with TCCs achieved an error rate of 2.59567387687% (2_2 labels).
# Spectral clustering with TCCs achieved an error rate of 38.7354409318% (9_2 labels).

# model = NMF(n_components=100, init='random', random_state=0)
# Spectral clustering with TCCs achieved an error rate of 2.59567387687% (2_2 labels).
# Spectral clustering with TCCs achieved an error rate of 38.6356073211% (9_2 labels).


W = model.fit_transform(D)
H = model.components_
print time.time() - a

D_2 = np.dot(W, H)
diff = D - D_2
print(diff)
print(D_2.shape)



# ---------------------------------------------
X_tsne = tSNE_pairwise(D)

print("Done with tsne")

# obtain labels via spectral clustering
def spectral(k,D):
    if D[1,1] < 1: D = 1-D # Convert distance to similarity matrix
    spectral = cluster.SpectralClustering(n_clusters=k,affinity='precomputed')
    spectral.fit(D)
    labels = spectral.labels_
    return labels

tcc_spectral_labels2 = spectral(2,D)
tcc_spectral_labels2_2 = spectral(2,D_2)
tcc_spectral_labels9 = spectral(9,D)
tcc_spectral_labels9_2 = spectral(9,D_2)
# tcc_spectral_labels47 = spectral(47,D)

# obtain labels via affinity propagation
def AffinityProp(D,pref,damp):
    aff= cluster.AffinityPropagation(affinity='precomputed',preference=pref,damping=damp, verbose=True)
    labels=aff.fit_predict(D)
    return labels

pref = -np.median(D.flatten())*np.ones(3005)
tcc_affinity_labels1 = AffinityProp(-D,pref,0.5)
tcc_affinity_labels2 = AffinityProp(-D,2*pref,0.7)


# gets max weight matching of a biparetite graph with row_label x column_label
# (weights are given by weight_matrix)
def get_max_wt_matching(row_label,column_label, weight_matrix):
    # Create a bipartite graph where each group has |unique labels| nodes 
    G = nx.complete_bipartite_graph(len(row_label), len(column_label))
    # Weight each edge by the weight in weight matrix.. 
    for u,v in G.edges(): G[u][v]["weight"]=weight_matrix[u,v-len(row_label)]
    # Perform weight matching using Kuhn Munkres
    H=nx.max_weight_matching(G)
    max_wt=0
    for u,v in H.items(): max_wt+=G[u][v]["weight"]/float(2)
    return max_wt

def compute_clustering_accuracy(label1, label2):
    uniq1,uniq2 = np.unique(label1),np.unique(label2)
    # Create two dictionaries. Each will store the indices of each label
    entries1,entries2 = {},{}
    for label in uniq1: entries1[label] = set(np.flatnonzero((label1==label)))
    for label in uniq2: entries2[label] = set(np.flatnonzero((label2==label)))
    # Create an intersection matrix which counts the number of entries that overlap for each label combination        
    W = np.zeros((len(uniq1),len(uniq2)))
    for i,j in itertools.product(range(len(uniq1)),range(len(uniq2))):
        W[i,j]=len(entries1[uniq1[i]].intersection(entries2[uniq2[j]]))
    # find the max weight matching
    match_val = get_max_wt_matching(uniq1,uniq2,W)
    # return the error rate
    return (1-match_val/float(len(label1)))*100

# c2 = compute_clustering_accuracy(tcc_spectral_labels2,labels2)
# print 'Spectral clustering with TCCs achieved an error rate of ' + str(c2) + '% (2 labels).'
# c9 = compute_clustering_accuracy(tcc_spectral_labels9,labels9)
# print 'Spectral clustering with TCCs achieved an error rate of ' + str(c9) + '% (9 labels).'
c2_2 = compute_clustering_accuracy(tcc_spectral_labels2_2,labels2)
print 'Spectral clustering with TCCs achieved an error rate of ' + str(c2_2) + '% (2_2 labels).'
c9_2 = compute_clustering_accuracy(tcc_spectral_labels9_2,labels9)
print 'Spectral clustering with TCCs achieved an error rate of ' + str(c9_2) + '% (9_2 labels).'



# N = len(D)
# M = len(D[0])
# K = 2

# P = np.random.rand(N,K)
# Q = np.random.rand(M,K)

# nP, nQ = matrix_factorization(D, P, Q, K)
# nR = np.dot(nP, nQ.T)
# print(nR)

