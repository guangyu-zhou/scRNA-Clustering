
import pickle 
import numpy as np
import time
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.manifold import TSNE
from sklearn import cluster,manifold
import matplotlib.pyplot as plt
# %matplotlib inline

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

def load_GE(path):
	dge = pd.read_csv(path + "processed_gene.csv", index_col = 0, sep = '\t')
	# Get the expression matrix: Cell * Genes
	X = dge.values.T
	cell_names = dge.columns.values

def load_21(path):
	with open(path + "Trapnell_TCC_21.dat.dat", 'rb') as f:
    	X=pickle.load(f, encoding='latin1')
    return X

def compute_dist():
	with open(path + 'Trapnell_TCC_pairwise_distance_dge.dat','wb') as outfile:
    	pickle.dump(D, outfile)

path = '/home/zgy_ucla_cs/Research/singleCell/TCC_old_pipeline/scRNA-Clustering/Trapnell_pipeline/'

with open(path + "Trapnell_TCC_pairwise_distance_21.dat", 'rb') as f:
    D=pickle.load(f, encoding='latin1')

with open(path + "Trapnell_TCC_pairwise_distance_31.dat", 'rb') as f:
    D=pickle.load(f, encoding='latin1')    

cluster_labels = np.loadtxt(path + 'Trapnells_data/Trapnell_labels.txt',dtype=str).astype(int)-1    
num_of_clusters=3
similarity_mat= D.max()-D
labels_spectral = spectral(num_of_clusters,similarity_mat)
print(NMI(cluster_labels, labels_spectral), ARI(cluster_labels, labels_spectral))  

# ===================== scVI =====================
# expression_train, expression_test, cluster_labels, c_test = train_test_split(X, X_type, random_state=0)


batch_size = 128
learning_rate = 0.001
epsilon = 0.01
latent_dimension = 10





tf.reset_default_graph()
expression = tf.placeholder(tf.float32, (None, X.shape[1]), name='x')
kl_scalar = tf.placeholder(tf.float32, (), name='kl_scalar')
optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate, epsilon=epsilon)
training_phase = tf.placeholder(tf.bool, (), name='training_phase')

# getting priors
log_library_size = np.log(np.sum(X, axis=1))
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
result = train_model(model, (X, X), sess, 250, batch_size=batch_size)

dic_full = {expression: X, training_phase:False}
latent = sess.run(model.z, feed_dict=dic_full)
# clustering_score = cluster_scores(latent, len(cell_types), cluster_labels)
clustering_score = cluster_scores(latent, np.max(cluster_labels), cluster_labels)
print("Silhouette", clustering_score[0], "\nAdjusted Rand Index", clustering_score[1], \
        "\nNormalized Mutual Information", clustering_score[2])