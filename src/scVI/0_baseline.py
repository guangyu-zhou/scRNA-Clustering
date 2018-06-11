import scVI
import tensorflow as tf
from benchmarking import *
from helper import *

import json
import numpy as np
import pickle
import pandas as pd
from sklearn import cluster,manifold
import matplotlib.pyplot as plt

import time
from sklearn.model_selection import train_test_split
from sklearn.manifold import TSNE

from sklearn.metrics import silhouette_score
from sklearn.cluster import KMeans
from sklearn.metrics import normalized_mutual_info_score as NMI
from sklearn.metrics import adjusted_rand_score as ARI

from sklearn.decomposition import TruncatedSVD
from sklearn.metrics.pairwise import pairwise_distances

def AffinityProp(D,pref,damp):
    aff= cluster.AffinityPropagation(affinity='precomputed',
                                     preference=pref,damping=damp, verbose=True)
    labels=aff.fit_predict(D)
    return labels

def spectral(k,D):
    spectral = cluster.SpectralClustering(n_clusters=k,affinity='precomputed')
    spectral.fit(D)
    labels = spectral.labels_
    return labels

# def spectral(k,D):
    # labels = cluster.spectral_clustering(D, n_clusters=k)
    # return labels    

def tSNE_pairwise(D):
    tsne = manifold.TSNE(n_components=2, random_state=213, metric='precomputed', n_iter=2000, verbose=1);
    X_tsne = tsne.fit_transform(D);
    return X_tsne


def load_Hong_DGE():
    dge = pd.read_csv("/home/zgy_ucla_cs/Research/DropSeq/data/dge_E31.csv", index_col = 0)
    X = dge.values.T
    print("X shape", X.shape)
    
    X_type = pd.read_csv("/home/zgy_ucla_cs/Research/DropSeq/data/cell.type_E31.csv", sep=",")
    cell_types, cluster_labels = np.unique(np.array(X_type)[:,1], return_inverse=True)
    print("Number of cluster labels",len(cluster_labels), len(cell_types))
    return X, cluster_labels
    num_of_clusters=len(cell_types)
    expression_train, expression_test, c_train, c_test = train_test_split(X, cluster_labels, random_state=0)
    # labels_spectral = spectral(expression_train, num_of_clusters)
    # print(NMI(c_train, labels_spectral), ARI(c_train, labels_spectral))
    return expression_train, expression_test, c_train, c_test



def load_Hong_21mer():
    path = '/home/zgy_ucla_cs/Research/singleCell/scRNA-Seq-TCC-prep/mat_21mer/'
    with open(path+"TCC_matrix.dat", 'rb') as f:
        # TCC=pickle.load(f, encoding='latin1')
        TCC=pickle.load(f)
        X = TCC.T
        print("non-zero ratio", 1- TCC.count_nonzero()*1.0/(TCC.shape[0]*TCC.shape[1]), TCC.nnz)
        print("X shape", X.shape)
    
    X_type = pd.read_csv("../../data/cell.type_E31.csv", sep=",")
    cell_types, cluster_labels = np.unique(np.array(X_type)[:,1], return_inverse=True)
    print("Number of cluster labels",len(cluster_labels), len(cell_types))
    
    return X, cluster_labels
    # with open(parameter["SAVE_DIR"]+"pwise_dist_L1.dat", 'rb') as f:
        # D_l1=pickle.load(f, encoding='latin1')
    # with open(parameter["SAVE_DIR"]+"nonzero_ec.dat", 'rb') as f:
        # nonzero_ec=pickle.load(f, encoding='latin1')

def load_Hong_31mer():
    path = '/home/zgy_ucla_cs/Research/singleCell/scRNA-Seq-TCC-prep/mat_31mer/'
    with open(path+"TCC_matrix.dat", 'rb') as f:
        # TCC=pickle.load(f, encoding='latin1')
        TCC=pickle.load(f)
        print("non-zero ratio", 1- TCC.count_nonzero()*1.0/(TCC.shape[0]*TCC.shape[1]), TCC.nnz)
        X = TCC.T
        print("X shape", X.shape)
    
    X_type = pd.read_csv("../../data/cell.type_E31.csv", sep=",")
    cell_types, cluster_labels = np.unique(np.array(X_type)[:,1], return_inverse=True)
    print("Number of cluster labels",len(cluster_labels), len(cell_types))
    
    return X, cluster_labels


def latent_dist(X, dim):
    pca = TruncatedSVD(n_components = dim)
    latent = pca.fit_transform(X)
    print("latent shape", latent.shape)
    num_of_threads = 8
    latent_dist = pairwise_distances(latent, metric='euclidean',n_jobs=num_of_threads)
    return latent_dist

''' Hong DGE '''
expression_train, expression_test, c_train, c_test = load_Hong_DGE()
path = '/home/zgy_ucla_cs/Research/DropSeq/data/'
with open(path+"dge_pwise_dist_JS.dat", 'rb') as f:
    D_l2 = pickle.load(f)#, encoding='latin1')

expression_train
num_of_clusters=8
similarity_mat=D_l2.max()-D_l2
labels_spectral = spectral(num_of_clusters,similarity_mat)
print(NMI(cluster_labels, labels_spectral), ARI(cluster_labels, labels_spectral))

''' Hong 21 '''
path = '/home/zgy_ucla_cs/Research/singleCell/scRNA-Seq-TCC-prep/mat_21mer/'
X, cluster_labels = load_Hong_21mer()
with open(path+"pwise_dist_JS.dat", 'rb') as f:
    D_l2=pickle.load(f, encoding='latin1')
num_of_clusters=8
similarity_mat=D_l2.max()-D_l2
labels_spectral = spectral(num_of_clusters,similarity_mat)
print(NMI(cluster_labels, labels_spectral), ARI(cluster_labels, labels_spectral))

pwise_dist_latent = latent_dist(X, 2)
similarity_mat=pwise_dist_latent.max()-pwise_dist_latent
labels_spectral = spectral(num_of_clusters,similarity_mat)
print(NMI(cluster_labels, labels_spectral), ARI(cluster_labels, labels_spectral))



''' Hong 31 '''
path = '/home/zgy_ucla_cs/Research/singleCell/scRNA-Seq-TCC-prep/mat_31mer/'
X, cluster_labels = load_Hong_31mer()
with open(path+"pwise_dist_JS.dat", 'rb') as f:
    D_js=pickle.load(f)#, encoding='latin1')

D = D_js
num_of_clusters=8
similarity_mat= D.max()-D
labels_spectral = spectral(num_of_clusters,similarity_mat)
print(NMI(cluster_labels, labels_spectral), ARI(cluster_labels, labels_spectral))

pwise_dist_latent = latent_dist(X, 10)
similarity_mat=pwise_dist_latent.max()-pwise_dist_latent
labels_spectral = spectral(num_of_clusters,similarity_mat)
print(NMI(cluster_labels, labels_spectral), ARI(cluster_labels, labels_spectral))


from scipy.stats.stats import pearsonr
from numpy import corrcoef

pca = TruncatedSVD(10)
latent = pca.fit_transform(X)
print(latent.shape)

X_d = X.todense()
Xn_d = []
Corn_d = []
for elem in [0,1,2,3,4,5,6,7]:
    ind_i = np.where(labels_spectral == elem)[0]
    print(ind_i.shape)
    X_elem = X_d[ind_i,:]
    Xn_d.append(X_elem)
    corn_elem = corrcoef(X_elem)
    Corn_d.append(corn_elem)
    print(elem, np.mean(corn_elem))





Xn_true = []
Corn_true = []
for elem in [0,1,2,3,4,5,6,7]:
    ind_i = np.where(cluster_labels == elem)[0]
    # labels_spectral
    print(ind_i.shape)
    X_elem = X[ind_i,:]
    Xn_true.append(X_elem)
    corn_elem = corrcoef(X_elem)
    Corn_true.append(corn_elem)
    print(elem, np.mean(corn_elem))

ind_i = np.where(labels_spectral == 0)[0]
labels_spectral[ind_i]=40

ind_i = np.where(labels_spectral == 1)[0]
labels_spectral[ind_i]=0

ind_i = np.where(labels_spectral == 2)[0]
labels_spectral[ind_i]=70

ind_i = np.where(labels_spectral == 4)[0]
labels_spectral[ind_i]=50

ind_i = np.where(labels_spectral == 5)[0]
labels_spectral[ind_i]=20

ind_i = np.where(labels_spectral == 6)[0]
labels_spectral[ind_i]=10

ind_i = np.where(labels_spectral == 7)[0]
labels_spectral[ind_i]=60

ind_i = np.where(labels_spectral == 10)[0]
labels_spectral[ind_i]=1

ind_i = np.where(labels_spectral == 20)[0]
labels_spectral[ind_i]=2

ind_i = np.where(labels_spectral == 30)[0]
labels_spectral[ind_i]=3

ind_i = np.where(labels_spectral == 40)[0]
labels_spectral[ind_i]=4

ind_i = np.where(labels_spectral == 50)[0]
labels_spectral[ind_i]=5

ind_i = np.where(labels_spectral == 60)[0]
labels_spectral[ind_i]=6

ind_i = np.where(labels_spectral == 70)[0]
labels_spectral[ind_i]=7

labels_spectral[ind_i]=7
labels_spectral_mapped = labels_spectral

from sklearn.metrics import confusion_matrix
y_true = [2, 0, 2, 2, 0, 1]
y_pred = [0, 0, 2, 2, 0, 2]
confusion_matrix(y_true, y_pred)

 0,4
 1,0
 2,7
 3,3
 4,5
 5,2
 6,1
 7,6
Pearson correlation, within Spectral label groups
(994,) 0.412490030222 [5]
(403,) 0.610216355785 [1?]
(256,) 0.41946210669 [8]
(317,) 0.389508953441 [4?]
(398,) 0.890224746369 [6]
(326,) 0.733754499755 [3]
(273,) 0.571539563932 [2]
(414,) 0.387209222224 [7]

Pearson correlation, within Ground truth groups
(918,) 0.649405603455 [1]
(200,) 0.463251764059 [2]
(332,) 0.706833689212 [3]
(77,) 0.408436931297 [4]
(970,) 0.326835582062 [5]
(393,) 0.895694268137 [6]
(413,) 0.407163835502 [7]
(78,) 0.443793744044 [8]
# spectral_clu = cluster.SpectralClustering(n_clusters=num_of_clusters)
# labels_spectral = spectral_clu.fit_predict(latent)





