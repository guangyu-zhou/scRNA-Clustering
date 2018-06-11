import matplotlib as mpl
mpl.use('Agg')
# import matplotlib.pyplot as plt
import numpy as np
import math
import pylab as plt
# import matplotlib
# matplotlib.use('Agg')
path = '/home/zgy_ucla_cs/data/dropSeq/'

with open(path + "E31_REP3_HHNKFBGX3_1.fastq") as infile:
	count = 0
	readname = ''
	barcode = ''
	bar_pool = set()
	bar_read = {}
	for line in infile:
		line = line.strip()
		count+=1
		if count % 1000000 == 0:
			print(count)
		if count % 4 == 1:
			readname = line
		if count % 4 == 2:
			barcode = line
		if count % 4 == 3:
			# print(barcode, readname)
			# bar_pool.add(barcode)
			if barcode not in bar_read:
				bar_read[barcode] = [readname]
			else:
				bar_read[barcode].append(readname)
		# if count > 1000000:
			# break


	print(len(bar_read))
	read_cluster_size = []
	for elem in bar_read.values():
		if len(elem) > 2:
			read_cluster_size.append(len(elem))
	# print read_cluster_size
	

	bins = np.linspace(math.ceil(min(read_cluster_size)), 
                   math.floor(max(read_cluster_size)),
                   20) # fixed number of bins

	
	plt.hist(read_cluster_size, bins = bins)  # arguments are passed to np.histogram
	plt.title("Histogram with 'auto' bins")
	plt.gca().set_yscale("log")
	plt.savefig("E31_REP3_HHNKFBGX3_1.png")
	# for k in bar_read:
		# print len(bar_read[k])
		# print(k, bar_read[k])