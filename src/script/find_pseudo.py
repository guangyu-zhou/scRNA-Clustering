import sys, re, os, argparse, datetime
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import collections

def load_pseudo_gene():
	pseudo_genes = set()
	with open("/home/zgy_ucla_cs/Research/DropSeq/data/annotation/pseudo_genes.gtf", 'r') as pseudo:
		for line in pseudo:
			l = line.strip().split(";")[0].split()[1]
			pseudo_genes.add(l.strip("\"")) 
	print("pseudo_genes loaded, size = ", len(pseudo_genes))

	return pseudo_genes


def load(outfile, pseudo_genes_freq = False):
	pseudo = load_pseudo_gene()
	total_mult_genes = set()
	total_pseduo_genes = set()
	cluster_cnt_multi = 0
	cluster_cnt_unique = 0
	pseudo_genes_count_arr = []
	# with open("/home/zgy_ucla_cs/Research/DropSeq/data/annotation/multi_annotations.txt", 'r') as ann:
	if outfile:
		print("writing to file: ", outfile)
		fout = open(outfile, "w")
		sys.stdout = fout
	with open("/home/zgy_ucla_cs/Research/DropSeq/data/annotation/majority.txt", 'r') as ann:
		for line in ann:
			max_count = 0
			is_unique = True
			l = line.strip().split()
			# print(l)
			valid_reads = int(l[1]) - int(l[4])
			# thre = valid_reads/2
			genes=l[-1].split(',')
			mult_genes = []
			for elem in genes:
				gene, sup = elem.split(':')
				# print sup, thre
				if int(sup)*2 >= valid_reads and gene!='NA':
					# if max_count < sup:
					# 	max_count = sup
					# 	# is_unique = True
					# elif max_count == sup:
					# 	is_unique = False

					mult_genes.append(gene)
					total_mult_genes.add(gene)
			# if len(mult_genes) == 0:
				# print l
			if len(mult_genes) < 2 :
				# if is_unique and len(mult_genes) >= 2:
					# print ("thre", thre)
					# print "is_unique"
					# print len(mult_genes), genes
				# cluster_cnt_unique+=1
				continue
			cluster_cnt_multi+=1
			# print l[0],
			mult_genes_pseudo = []
			for g in mult_genes:
				if g in pseudo:
					mult_genes_pseudo.append(g)
					total_pseduo_genes.add(g)
				# print (g, g in pseudo),
			
			# if len(mult_genes_pseudo) > 0:
				# print l[0], (len(mult_genes), len(mult_genes_pseudo), mult_genes_pseudo)
			if pseudo_genes_freq:
				pseudo_genes_count_arr.append(len(mult_genes_pseudo))
			# break

	print("total_mult_genes", "total_pseduo_genes", len(total_mult_genes), len(total_pseduo_genes))
	print("#clusters with mult_genes", cluster_cnt_multi)

	sys.stdout = sys.__stdout__

	# print(pseudo_genes_count_arr)
	return pseudo_genes_count_arr

def main(parser):
	options = parser.parse_args()
	outfile = options.outfile
	
	pseudo_genes_count_arr = load(outfile, pseudo_genes_freq = True)
	counter=collections.Counter(pseudo_genes_count_arr)
	print counter
	# if len(pseudo_genes_count_arr):
		# plt.hist(pseudo_genes_count_arr)
		# plt.savefig("tmp.jpg")


if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog="find_pseudo.py")
	parser.add_argument("-o", "--outfile", dest="outfile", type=str, help="output pseudo genes stats filename", required = False)
	main(parser)





