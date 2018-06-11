import itertools
import matplotlib as mpl
mpl.use('Agg')
# import matplotlib.pyplot as plt
import numpy as np
import math
import pylab as plt
import pickle

# import matplotlib
# matplotlib.use('Agg')

def split_fastq(path, barcode_arr):
	
	with open(path + "E31_REP3.fastq", 'r') as f:
		count = 0
		for line1,line2, line3, line4 in itertools.izip_longest(*[f]*4):
			if barcode_arr[count]!='-':
				print(count, barcode_arr[count])
				# fout = open(path + '/cells/' + barcode_arr[count][0] + '.fastq', 'a+')
				# fout.write(line1)
				# fout.write(line2)
				# fout.write(line3)
				# fout.write(line4)

				fout_umi = open(path + '/cells/' + barcode_arr[count][0] + '.umi', 'a+')
				fout_umi.write(barcode_arr[count][1] + '\n')
				
			count+=1
			# break

def load_barcode_arr(path):
	# with open(path + "E31_REP3.barcode", 'r') as infile:
	# 	count = 0
	# 	readname = ''
	# 	barcode = ''
	# 	barcode_arr = []
	# 	cell_read = {}
	# 	for line in infile:
	# 		readname,barcode = line.strip().split('\t')
	# 		cell_barcode = barcode[:12]
	# 		umi = barcode[12:]
	# 		barcode_arr.append([cell_barcode, umi])
	# 		if cell_barcode not in cell_read:
	# 			cell_read[cell_barcode] = [readname]
	# 		else:
	# 			cell_read[cell_barcode].append(readname)
			
	# 		# print(readname,cell_barcode)
	# 	print("Done with reading, start filtering")
	# 	for i in range(len(barcode_arr)):
	# 		readCt = len(cell_read[barcode_arr[i][0]]) 
	# 		# if readCt < 2 or readCt > 100:
	# 		if readCt != 50:
	# 			barcode_arr[i] = '-'
	# 	print(len(barcode_arr))
	# 	pickle.dump( barcode_arr, open( path + "barcode_arr_50.p", "wb" ) )
	barcode_arr = pickle.load(open( path + "barcode_arr_50.p", "rb" ))
	return barcode_arr
		

# with open(path + "E31_REP3.barcode", 'r') as infile:
# 		count = 0
# 		line_cnt = 0
# 		readname = ''
# 		barcode = ''
# 		bar_pool = set()
# 		cell_read = {}
# 		for line in infile:
# 			readname,barcode = line.strip().split('\t')
# 			cell_barcode = barcode[:12]
# 			count+=1

# 			if cell_barcode not in cell_read:
# 				cell_read[cell_barcode] = [readname]
# 			else:
# 				cell_read[cell_barcode].append(readname)
# 			print(readname,cell_barcode)

# 			if count > 1000:
# 				break


# 		print(len(cell_read))
# read_cluster_size = []
# for elem in cell_read.values():
# 	if len(elem) > 2:
# 		read_cluster_size.append(len(elem))
# print "Cells with more than 2 elems",len(read_cluster_size)

# np.save(path + 'E31_REP3_cell_group_size', np.array(read_cluster_size))		

path = '/home/zgy_ucla_cs/Research/DropSeq/data/reads_barcode/'

barcode_arr = load_barcode_arr(path)
print("Done loading")
# for elem in barcode_arr:
# 	if elem!='-':
# 		print elem
# barcode_arr = ['aaa']
split_fastq(path, barcode_arr)
	