import numpy as np
import math
import pickle

# function that return a index file: key readname, value barcode (full length)
def build_umi_index(path):
	with open(path + "E31_REP2_HHN7NBGX3_1.fastq", "r") as infile:
		count = 0
		readname = ''
		barcode = ''
		read2bar = {}
		for line in infile:
			line = line.strip()
			count+=1
			if count % 1000000 == 0:
				print(count)
			elif count % 4 == 1:
				readname = line.split()[0]
			elif count % 4 == 2:
				barcode = line
			elif count % 4 == 3:
				# print(barcode, readname)
				read2bar[readname]=barcode
			# if count > 100:
				# break
	
	print("Dumping file") # use stderr to print out instant msg
	pickle.dump( read2bar, open( path + "E31_REP2_index.p", "wb" ) )
	return read2bar

path = '/home/zgy_ucla_cs/Research/DropSeq/data/Reps/'
read2bar = build_umi_index(path)
with open(path + "E31_REP2_HHN7NBGX3_reads_filtered.fastq") as infile:
	fout = open(path + "E31_REP2.umi", 'w')
	print("Starting to match umi")
	for line in infile:
		line = line.strip()
		# line = '@NS500551:319:HKT73BGX2:1:11101:17957:1049'
		umi = read2bar[line]
		fout.write(line + '\t' +umi + '\n')
		# print(line, umi)
		# break






