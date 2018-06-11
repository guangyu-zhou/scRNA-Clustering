import itertools
import os 

path = '/home/zgy_ucla_cs/Research/DropSeq/data/reads_barcode/cells/'

flnames=sorted([x for x in os.listdir(path) if x.endswith('.fastq')])
for flname in flnames:
	with open(path + flname, 'r') as f:
		cnt = 0
		flname = flname.split('.')[0]
		print flname
		fout_umi = open('/home/zgy_ucla_cs/Research/DropSeq/data/reads_barcode/cell_num_all/' + flname + '.fastq', 'w')
		for line1,line2, line3, line4 in itertools.izip_longest(*[f]*4):
			cnt+=1
			line1_new = '@'+flname + '.' + str(cnt)
			fout_umi.write(line1_new + '\n')
			fout_umi.write(line2)
			fout_umi.write(line3)
			fout_umi.write(line4)