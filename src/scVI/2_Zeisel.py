import scVI
import tensorflow as tf
from benchmarking import *
from helper import *

import numpy as np
import time
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.manifold import TSNE
from sklearn import cluster,manifold
import matplotlib.pyplot as plt
# %matplotlib inline
import pickle 

from sklearn.metrics import silhouette_score
from sklearn.cluster import KMeans
from sklearn.metrics import normalized_mutual_info_score as NMI
from sklearn.metrics import adjusted_rand_score as ARI

from sklearn.decomposition import TruncatedSVD
from sklearn.metrics.pairwise import pairwise_distances

from scipy.stats.stats import pearsonr
from numpy import corrcoef

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


def load_Zeisel_kmer31():
    path0 = "/home/zgy_ucla_cs/Research/singleCell/TCC_old_pipeline/scRNA-Clustering/Zeisel_pipeline/mat_31/"
    path = '/home/zgy_ucla_cs/Research/singleCell/scRNA-Seq-TCC-prep/Zeisel/'
    with open(path0+"TCC_matrix.dat", 'rb') as f:
        T=pickle.load(f)#, encoding='latin1')
        X = T.T
    with open(path0+"pwise_dist_L1.dat", 'rb') as f:
        D_l1=pickle.load(f)#, encoding='latin1')
    labels9 = np.loadtxt(path + 'Zeisels_labels9.txt',dtype=str).astype(int)-1    
    return X, labels9

def load_Zeisel_GE():
    path = '/home/zgy_ucla_cs/Research/singleCell/scRNA-Seq-TCC-prep/Zeisel/'
    dge = pd.read_csv(path + "GSE60361_C1-3005-Expression.txt", sep = '\t' ,index_col = 0)
    print(dge.shape)
    X = dge.values.T
    labels9 = np.loadtxt(path + 'Zeisels_labels9.txt',dtype=str).astype(int)-1
    return X, labels9
    # labels9 = labels9[:3004]
    # Zeisel's neurons (labeled as 0) and non-neurons (labeled as 1) obtained from labels9
    # labels2 = np.copy(labels9)
    # for i in [0,1,2]: labels2[labels2 == i] = 0
    # labels2[labels2 != 0] = 1
def load_Zeisel_GE2():
    path2 = '/home/zgy_ucla_cs/Research/singleCell/scRNA-Seq-TCC-prep/Zeisel2/'
    X2 = pd.read_csv(path2 + "expression_mRNA_17-Aug-2014.txt", sep="\t", low_memory=False).T
    clusters = np.array(X2[7], dtype=str)[2:]
    cell_types, labels = np.unique(clusters, return_inverse=True)
    gene_names = np.array(X.iloc[0], dtype=str)[10:]
    X = X.loc[:, 10:]
    X = X.drop(X.index[0])
    expression_data = np.array(X, dtype=np.int)[1:]

    # keep the most variable genes according to the Biscuit ICML paper
    selected = np.std(expression_data, axis=0).argsort()[-558:][::-1]
    expression_data = expression_data[:, selected]
    gene_names = gene_names[selected].astype(str)

def check_cluster_cor(cluster_labels):
    Xn_true = []
    Corn_true = []
    for elem in np.unique(cluster_labels):
        ind_i = np.where(cluster_labels == elem)[0]
        # labels_spectral
        print(ind_i.shape)
        X_elem = X[ind_i,:]
        Xn_true.append(X_elem)
        corn_elem = corrcoef(X_elem)
        Corn_true.append(corn_elem)
        print(elem, np.mean(corn_elem))
def check_cluster_cor_spectral():
    Xn_d = []
    Corn_d = []
    for elem in [0,1,2,3,4,5,6,7,8]:
        ind_i = np.where(labels_spectral == elem)[0]
        print(ind_i.shape)
        X_elem = X[ind_i,:]
        Xn_d.append(X_elem)
        corn_elem = corrcoef(X_elem)
        Corn_d.append(corn_elem)
        print(elem, np.mean(corn_elem))


#train test split for log-likelihood scores
expression_train, expression_test, c_train, c_test = train_test_split(expression_data, labels, random_state=0)
    # pickle.dump(expression_data, open(path + 'dge.dat', "wb"))

X, cluster_label = load_Zeisel_GE()

X, X_type = load_Zeisel_kmer31()
X = X.toarray()


# ============================== spectral ==============================

# Load Z TCC 31
path = '/home/zgy_ucla_cs/Research/singleCell/TCC_old_pipeline/scRNA-Clustering/Zeisel_pipeline/mat_31/'
with open(path+"pwise_dist_L1.dat", 'rb') as f:
    D=pickle.load(f)#, encoding='latin1')

path = '/home/zgy_ucla_cs/Research/singleCell/scRNA-Seq-TCC-prep/Zeisel/'
with open(path+"pwise_dist_l1.dat", 'rb') as f:
    D_l1=pickle.load(f)#, encoding='latin1')

D = D_l1
num_of_clusters=9
similarity_mat= D.max()-D
labels_spectral = spectral(num_of_clusters,similarity_mat)
print(NMI(cluster_labels, labels_spectral), ARI(cluster_labels, labels_spectral))

# pwise_dist_latent = latent_dist(X, 10)
# similarity_mat=pwise_dist_latent.max()-pwise_dist_latent
# labels_spectral = spectral(num_of_clusters,similarity_mat)
# print(NMI(cluster_labels, labels_spectral), ARI(cluster_labels, labels_spectral))






# ===================== scVI =====================
# expression_train, expression_test, c_train, c_test = train_test_split(X, X_type, random_state=0)
expression_train, expression_test, c_train, c_test = train_test_split(expression_data, X_type, random_state=0)

log_library_size = np.log(np.sum(expression_train, axis=1))
mean, var = np.mean(log_library_size), np.var(log_library_size)

batch_size = 128
learning_rate = 0.001
epsilon = 0.01
latent_dimension = 10





tf.reset_default_graph()
expression = tf.placeholder(tf.float32, (None, expression_train.shape[1]), name='x')
kl_scalar = tf.placeholder(tf.float32, (), name='kl_scalar')
optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate, epsilon=epsilon)
training_phase = tf.placeholder(tf.bool, (), name='training_phase')

# getting priors
log_library_size = np.log(np.sum(expression_train, axis=1))
mean, var = np.mean(log_library_size), np.var(log_library_size)

# loading data
model = scVI.scVIModel(expression=expression, kl_scale=kl_scalar, \
                         optimize_algo=optimizer, phase=training_phase, \
                          library_size_mean=mean, library_size_var=var, n_latent=latent_dimension)

#starting computing session
sess = tf.Session()

 # Initialize the graph and fit the training set
# this takes less than a minute on a Tesla K80
sess.run(tf.global_variables_initializer())
result = train_model(model, (expression_train, expression_test), sess, 250, batch_size=batch_size)

dic_full = {expression: expression_train, training_phase:False}
latent = sess.run(model.z, feed_dict=dic_full)
# clustering_score = cluster_scores(latent, len(cell_types), c_train)
clustering_score = cluster_scores(latent, np.max(c_train), c_train)
print("Silhouette", clustering_score[0], "\nAdjusted Rand Index", clustering_score[1], \
        "\nNormalized Mutual Information", clustering_score[2])