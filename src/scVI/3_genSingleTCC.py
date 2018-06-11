import json
import numpy as np
import pickle
with open("/home/zgy_ucla_cs/Research/singleCell/scRNA-Seq-TCC-prep/example_dataset/config_31.json") as json_file:
    parameter = json.load(json_file)
    
with open(parameter["SAVE_DIR"]+"TCC_matrix.dat", 'rb') as f:
    T=pickle.load(f, encoding='latin1')
with open(parameter["SAVE_DIR"]+"pwise_dist_L1.dat", 'rb') as f:
    D_l1=pickle.load(f, encoding='latin1')
with open(parameter["SAVE_DIR"]+"nonzero_ec.dat", 'rb') as f:
    nonzero_ec=pickle.load(f, encoding='latin1')

ecfile_dir = parameter["kallisto"]["TCC_output"]+'matrix.ec'
eclist=np.loadtxt(ecfile_dir,dtype=str)    
print("Done loading")

from sklearn.preprocessing import normalize

TCC=T.T
T_norm = normalize(T, norm='l1', axis=0) 
print("Done normalize")
T_normT = T_norm.transpose()
    
NUM_OF_CELLS=np.shape(T)[1]
print("NUM_OF_CELLS =", NUM_OF_CELLS)
print("NUM_OF_nonzero_EC =", np.shape(T)[0])    

EC_dict = {}
for i in range(np.shape(eclist)[0]):
    EC_dict[i] = [int(x) for x in eclist[i,1].split(',') if x != ''] # fix by add if statement
    
union=set()
for i in nonzero_ec:
    new = [tx for tx in EC_dict[i] if tx not in union] # filter out previously seen transcripts
#     if new != []:
#         print(new)
    union.update(new) 

union_list=list(union) #union of all transctipt ids seen in nonzero eq.classes
NUM_OF_TX_inTCC = len(union)
print("NUM_OF_DISTINCT_Transcripts =", NUM_OF_TX_inTCC) #number of distinct transcripts in nonzero eq. classes

TCC_singleT = T_normT[:,:NUM_OF_TX_inTCC]
print("TCC_singleT: ", TCC_singleT.shape)
print("non-zero ratio", 1- TCC_singleT.count_nonzero()*1.0/(TCC_singleT.shape[0]*TCC_singleT.shape[1]), TCC_singleT.nnz)

with open(parameter["SAVE_DIR"]+"TCC_singleT.dat", 'wb') as f:
    pickle.dump(TCC_singleT,f)



