# awk -F" " '{print $1}' E31_REP1_HKT73BGX2_1.fastq | grep '^@' > E31_REP1_HKT73BGX2_all_reads.fastq
# diff --new-line-format="" --unchanged-line-format="" <(sort E31_REP1_HKT73BGX2_all_reads.fastq) <(sort E31_REP1_HKT73BGX2_reads_filtered.fastq) > E31_REP1_HKT73BGX2_removed.fastq

import matplotlib as mpl
mpl.use('Agg')
# import matplotlib.pyplot as plt
import numpy as np
import math
import pylab as plt

def process_individual():
	path = '/home/zgy_ucla_cs/data/dropSeq/'

	low_qual_reads = set()

	with open(path + "E31_REP1_HKT73BGX2_removed.fastq") as infile0:
		for line in infile0:
			line = line.strip()
			low_qual_reads.add(line)
			# break

	print ("done with low_qual_reads reading")
	# print(len(low_qual_reads), low_qual_reads)

	bar_read = {}
	with open(path + "E31_REP1_HKT73BGX2_1.fastq") as infile:
		count = 0
		readname = ''
		barcode = ''
		for line in infile:
			line = line.strip()
			count+=1
			if count % 1000000 == 0:
				print(count)
			if count % 4 == 1:
				readname = line.split()[0]
				# print(readname)
			if count % 4 == 2:
				barcode = line
			if count % 4 == 3 and readname not in low_qual_reads:
			
			# if count % 4 == 3:
				# if readname in low_qual_reads:
					# print readname
				if barcode not in bar_read:
					bar_read[barcode] = [readname]
				else:
					bar_read[barcode].append(readname)
			# if count > 100000:
				# break

	print("Bar size", len(bar_read))
	f1 = open(path + "E31_REP1_HKT73BGX2_cluster_1.fastq", 'w+')
	f_other = open(path + "E31_REP1_HKT73BGX2_cluster_other.fastq", 'w+')
	for key in bar_read.keys():
		list1 = bar_read[key]
		str1 = ' '.join(list1)
		if len(list1) == 1:
			f1.write(key + ' ' + str1 + '\n')
		else:
			f_other.write(key + ' ' + str1 + '\n')

	print("start plotting histogram")
	read_cluster_size = []
	for elem in bar_read.values():
		if len(elem) > 1:
			read_cluster_size.append(len(elem))
	# print read_cluster_size


	bins = np.linspace(math.ceil(min(read_cluster_size)), 
	               math.floor(max(read_cluster_size)),
	               100) # fixed number of bins


	plt.hist(read_cluster_size, bins = bins)  # arguments are passed to np.histogram
	plt.title("Histogram with 100 bins, log scale")
	plt.gca().set_yscale("log")
	plt.savefig(path + "E31_REP1_HKT73BGX2_cluster_size.png")	

def process_all():
	path = '/home/zgy_ucla_cs/data/dropSeq/'
	low_qual_reads = set()
	rep_names = ['E31_REP1_HKT73BGX2', 'E31_REP2_HHN7NBGX3', 'E31_REP3_HHNKFBGX3']
	
	for rep in rep_names:
		with open(path + rep+ "_removed.fastq") as infile0:
			for line in infile0:
				line = line.strip()
				low_qual_reads.add(line)
				# break

	print ("done with low_qual_reads reading", len(low_qual_reads))

	bar_read = {}

	for rep in rep_names:
		with open(path + rep + "_1.fastq") as infile:
			print("processing", rep)
			count = 0
			readname = ''
			barcode = ''
			for line in infile:
				line = line.strip()
				count+=1
				if count % 1000000 == 0:
					print(count)
				if count % 4 == 1:
					readname = line.split()[0]
					# print(readname)
				if count % 4 == 2:
					barcode = line
				if count % 4 == 3 and readname not in low_qual_reads:
				
				# if count % 4 == 3:
					# if readname in low_qual_reads:
						# print readname
					if barcode not in bar_read:
						bar_read[barcode] = [readname]
					else:
						bar_read[barcode].append(readname)
				# if count > 100000:
					# break

		print("Bar size", len(bar_read))
		
	f1 = open(path + "All_cluster_1.fastq", 'w+')
	f_other = open(path + "All_cluster_other.fastq", 'w+')
	for key in bar_read.keys():
		list1 = bar_read[key]
		str1 = ' '.join(list1)
		if len(list1) == 1:
			f1.write(key + ' ' + str1 + '\n')
		else:
			f_other.write(key + ' ' + str1 + '\n')

	print("start plotting histogram")
	read_cluster_size = []
	for elem in bar_read.values():
		if len(elem) > 0:
			read_cluster_size.append(len(elem))
	# print read_cluster_size
	np.save(path + 'All_cluster_size', np.array(read_cluster_size))

	bins = np.linspace(math.ceil(min(read_cluster_size)), 
	               math.floor(max(read_cluster_size)),
	               500) # fixed number of bins


	plt.hist(read_cluster_size, bins = bins)  # arguments are passed to np.histogram
	plt.title("Cluster size distritbuion, log scale")
	plt.gca().set_yscale("log")
	plt.savefig(path + "All_cluster_size.png")	
process_all()

