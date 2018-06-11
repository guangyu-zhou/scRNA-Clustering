import scVI
import tensorflow as tf
from benchmarking import *
from helper import *

import numpy as np
import time
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
# %matplotlib inline

X_type = pd.read_csv("../../data/cell.type_E31.csv", sep=",")
cell_types, cluster_labels = np.unique(np.array(X_type)[:,1], return_inverse=True)
# le = preprocessing.LabelEncoder()
# le.fit(cell_types)
# cluster_labels = le.transform(cell_types)
print(len(cluster_labels))

# dge = np.load("../../data/dge_E31.dat")
dge = pd.read_csv("../../data/dge_E31.csv", index_col = 0)
# Get the expression matrix: Cell * Genes

X = dge.values.T
cell_names = dge.columns.values
cell_names2 = np.array(X_type)[:,0]

# Make sure the cell name matchs
for i in range(len(cell_names)):
	# print(cell_names[i], cell_names2[i])
	if cell_names[i] not in cell_names2[i]:
		print("Not match", cell_names[i], cell_names2)[i]

gene_names = dge.index.tolist()
print(dge.shape)

# split dataset into train and test
expression_train, expression_test, c_train, c_test = train_test_split(X, cluster_labels, random_state=0)
print(cell_types)
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
clustering_score = cluster_scores(latent, 9, c_train)
print ("Silhouette", clustering_score[0], "Normalized Mutual Information", clustering_score[2],"Adjusted Rand Index", clustering_score[1])
        "\n)