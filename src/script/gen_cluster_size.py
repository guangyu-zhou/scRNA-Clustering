import matplotlib as mpl
mpl.use('Agg')
# import matplotlib.pyplot as plt
import numpy as np
import math
import pylab as plt
from matplotlib.ticker import ScalarFormatter


# path = '/home/zgy_ucla_cs/data/dropSeq/'
path = '/home/zgy_ucla_cs/Research/DropSeq/data/reads_barcode/'

read_cluster_size = []
# with open(path + "E31_REP1_HKT73BGX2_cluster_other.fastq") as infile:
# with open(path + "E31_REP2_HHN7NBGX3_cluster_other.fastq") as infile: 
# with open(path + "E31_REP3_HHNKFBGX3_cluster_other.fastq") as infile: 
	# l = 10000
# 	for line in infile:
# 		# l-=1
# 		# print line.count('@')
# 		size = line.count('@')
# 		if size <= 100:
# 			read_cluster_size.append(size)
# np.save(path + 'arr2_100', np.array(read_cluster_size))

# np.save(path + 'All_cluster_size', np.array(read_cluster_size))
read_cluster_size = np.load(path + 'E31_REP3_cell_group_size.npy')
read_cluster_size = read_cluster_size[np.where( read_cluster_size == 50)]
# read_cluster_size = read_cluster_size[np.where( read_cluster_size <= 60)]
print("total clusters:", len(read_cluster_size))
# print(min(read_cluster_size), max(read_cluster_size))

# bins = np.linspace(math.ceil(min(read_cluster_size)), 
#                math.floor(max(read_cluster_size)),
#                200) # fixed number of bins
# plt.hist(read_cluster_size, bins = bins)  # arguments are passed to np.histogram
# plt.title("Cluster size distritbuion, y log scale")
# plt.xlabel("reads counts")
# # plt.gca().set_xscale("log")

# # tk = np.linspace(math.ceil(min(read_cluster_size)), 
#                # math.floor(max(read_cluster_size)),
#                # 20)
# # plt.gca().set_xticks(tk)
# plt.gca().set_yscale("log")
# plt.ylabel("number of cells")

# # for axis in [plt.gca().xaxis, plt.gca().yaxis]:
#     # axis.set_major_formatter(ScalarFormatter())
# plt.savefig(path + "E31_REP3_cell_group_size_2_more_100_less.png")


