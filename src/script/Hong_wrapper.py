# Runs the full Zeisel pipeline. 
# Start from a folder of downloaded SRR files. Zeisel's 3005 mouse brain cell dataset can be found at http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60361.
# The pipeline runs kallisto to get the TCCs for each cell, ultimately outputting the TCC matrix (3005-by-#Eq classes).
# The pipeline also uses the TCC matrix to compute the pairwise distances between cells, resulting in a 3005-by-3005 distance matrix.

import os
import getopt
import sys

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:n:k:t:",["idir=","njobs=","hacked-kallisto-path=","reference-transcriptome="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python Zeisel_wrapper.py -i input_SRA_dir -k path-to-hacked-kallisto -t path-to-human-reference-transcriptome [-n number-of-processes-to-use]')
    sys.exit(1)
    
SRA_dir=''
num_proc=1
kallipso_path=''
ref_transcriptome=''

for opt,arg in opts:
    if opt in ("-i", "--idir"):
        SRA_dir=arg
    elif opt in ("-n","--njobs"):
        num_proc=int(arg)
    elif opt in ("-k","--hacked-kallisto-path"):
        kallipso_path=arg
    elif opt in ("-t","--reference-transcriptome"):
        ref_transcriptome=arg
    
        
if (not SRA_dir) or (not kallipso_path) or (not ref_transcriptome):
    print ('usage is : \n python Zeisel_wrapper.py -i input_SRA_dir -k path-to-hacked-kallisto -t path-to-mouse-reference-transcriptome [-n number-of-processes-to-use]')
    sys.exit(1)

'''
# GY comment out
print('Extracting reads from SRAs...')
os.system('mkdir -p ./reads_with_UMIs/')
os.system('rm -f ./reads_with_UMIs/*')
os.system('python process_SRA.py -i '+SRA_dir+' -o ./reads_with_UMIs/ -n '+str(num_proc))
'''

# read_dir_base='/home/zgy_ucla_cs/Research/DropSeq/data/reads_barcode/cells/'
read_dir_base='/home/zgy_ucla_cs/Research/DropSeq/data/reads_barcode/cell_num_all/'
# read_dir_base = '/home/zgy_ucla_cs/Research/singleCell/TCC/scRNA-Clustering/Zeisel_pipeline/read_UMI_small/'
index_path=ref_transcriptome


# print('Generating the Kallisto index (with hacked kallisto)...')
# os.system('mkdir -p ./kallisto_index')
# os.system('rm -f ./kallisto_index/*')
# os.system(kallipso_path+' index -i '+index_path+' '+ref_transcriptome)
# metadata_cmd=kallipso_path+' metadata '+index_path
# os.system(metadata_cmd)
# num_ec = sum(1 for line in open('./kallisto_index/Zeisel_index.idx_ecmap.txt'))
# print(num_ec)


# 3-22
# print('Generating TCC (with hacked kallisto)...')
# TCC_base_dir='./transcript_compatibility_counts_subsample'
# for index in range(1):
#     print('Running hacked kallisto on '+sampling_rates[index]+' fraction of reads...')
#     TCC_dir=TCC_base_dir+sampling_suffix[index]+"/"
#     read_dir_to_pass=read_dir_base+sampling_suffix[index]+"/"
#     os.system('mkdir -p '+TCC_dir)
#     os.system('rm -f '+TCC_dir+'*')
#     os.system('python get_pseudoalignments.py -i '+read_dir_to_pass+' -o '+TCC_dir+' -k '+kallipso_path+ ' -t '+ index_path +' -n '+ str(num_proc))

# print('Generating TCC distribution...')
# TCC_dist_base_flname='./Zeisel_TCC_distribution_subsample'
# TCC_base_flname='./Zeisel_TCC_subsample'
# for index in range(1):
#     print('Getting the TCC dist for '+sampling_rates[index]+' fraction of reads...')
#     TCC_dir=TCC_base_dir+sampling_suffix[index]+"/"
#     TCC_dist_flname=TCC_dist_base_flname+sampling_suffix[index]+".dat"
#     TCC_flname=TCC_base_flname+sampling_suffix[index]+".dat"
#     os.system('python get_tcc_dist.py -i '+TCC_dir+' -m '+str(num_ec)+' -t '+TCC_flname+' -d '+ TCC_dist_flname)


# print('Generating pairwise distances...')
# TCC_distance_base_flname='Zeisel_TCC_pairwise_JS_distance_subsample'
# for index in range(1):
#     TCC_dist_flname=TCC_dist_base_flname+sampling_suffix[index]+".dat"
#     TCC_distance_flname=TCC_distance_base_flname+sampling_suffix[index]+".dat"
#     print('Getting  pairwise distances for '+sampling_rates[index]+' fraction of reads...')
#     os.system('python get_pairwise_distances.py '+TCC_dist_flname+' '+TCC_distance_flname+' '+str(num_proc))

# Done....

# ================================================================
# ''' Running Kallisto quant and writing out pseudo-bams '''
print('Running Kallisto quant and writing out pseudo-bams...')
quant_dir_base='kallisto_quant_50'
pbam_dir_base='pbam_50'
tmp_dir='./temp_dir/'
os.system('mkdir -p '+tmp_dir)
print('Running kallisto pbam and quant')
pbam_dir= pbam_dir_base+"/"
quant_dir= quant_dir_base+"/"
read_dir_to_pass=read_dir_base+"/"
os.system('python /home/zgy_ucla_cs/Research/singleCell/TCC/scRNA-Clustering/Zeisel_pipeline/run_kallisto_and_get_pbam.py -i '+read_dir_to_pass+' -o '+quant_dir+' -p '+pbam_dir+' -t '+tmp_dir+' -r '+index_path +' -n '+str(num_proc))

# ================================================================    

pbam_dir_base='pbam_50'
TCC_UMI_tmp_base_dir='./Hong_TCC_UMI_tmp_50'
UMI_base_dir=read_dir_base+"/"
TCC_UMI_tmp_dir=TCC_UMI_tmp_base_dir+"/"
pbam_dir= pbam_dir_base+"/"
os.system('python /home/zgy_ucla_cs/Research/singleCell/TCC/scRNA-Clustering/Zeisel_pipeline/process_pbam.py -i '+pbam_dir+' -o '+TCC_UMI_tmp_dir+' -u '+UMI_base_dir+' -n '+str(num_proc))


TCC_UMI_file_base='./Hong_TCC_UMI_50'
TCC_UMI_distribution_file_base='./Hong_TCC_UMI_distribution_50'
# for index in range(1):
    # print('Getting TCC UMI matrix for '+sampling_rates[index]+' fraction of reads...')
TCC_UMI_tmp_dir=TCC_UMI_tmp_base_dir+"/"
TCC_UMI_file=TCC_UMI_file_base+'.dat'
TCC_UMI_distribution_file=TCC_UMI_distribution_file_base+'.dat'
os.system('python /home/zgy_ucla_cs/Research/singleCell/TCC/scRNA-Clustering/Zeisel_pipeline/get_UMI_matrices.py -i '+TCC_UMI_tmp_dir+' -t '+TCC_UMI_file+' -d '+TCC_UMI_distribution_file)
    
print('Generating pairwise distances...')
TCC_UMI_distance_file_base='./Hong_TCC_UMI_pairwise_JS_50'
# for index in range(1):
    # print('Getting pairwise distance of the TCC UMI matrix for '+sampling_rates[index]+' fraction of reads...')
TCC_UMI_distance_file=TCC_UMI_distance_file_base+'.dat'
TCC_UMI_distribution_file=TCC_UMI_distribution_file_base+'.dat'
os.system('python /home/zgy_ucla_cs/Research/singleCell/TCC/scRNA-Clustering/Zeisel_pipeline/get_pairwise_distances.py '+TCC_UMI_distribution_file+' '+TCC_UMI_distance_file+' '+str(num_proc))

        # print('Getting bowtie indices...')
        # bowtie_index_dir='./bowtie_index/'
        # os.system('mkdir -p '+bowtie_index_dir)
        # bowtie_index_path=bowtie_index_dir+'Zeisel_index.all'
        # os.system('bowtie-build --offrate=1 '+ref_transcriptome+' '+bowtie_index_path)

        # print('Running bowtie... May take days!!!')
        # read_dir_to_pass=read_dir_base+sampling_suffix[0]+"/"
        # bowtie_dir_base='./Zeisel_Bowtie_subsample'
        # bowtie_dir100='./Zeisel_Bowtie_subsample'+sampling_suffix[0]+"/"
        # os.system('python run_bowtie.py -i '+read_dir_to_pass+' -o '+bowtie_dir100+' -r '+bowtie_index_path+' -n '+str(num_proc))

        # print('Sampling bam files ...')
        # os.system('python subsample_bams.py -s '+read_dir_base+' -b ' +bowtie_dir100 +' -o '+bowtie_dir_base+' -n '+str(num_proc))

        # print('Counting UMIs...')
        # os.system('python UMI_counting.py  -n '+str(num_proc))
        # os.system('python UMI_counting_100.py  -n '+str(num_proc))

        # UMI_dir_base='./Zeisel_UMI_counts_subsample'
        # UMI_distribution_file_base='./Zeisel_UMI_distribution_subsample'
        # UMI_gene_file_base='./Zeisel_UMI_gene_distribution_subsample'
        # UMI_file_base='./Zeisel_UMI_subsample'
        # print('Getting UMI matrices....')
        # for index in range(1):
        #     UMI_file=UMI_file_base+sampling_suffix[index]+'.dat'
        #     UMI_distribution_file=UMI_distribution_file_base+sampling_suffix[index]+'.dat'
        #     UMI_gene_file=UMI_gene_file_base+sampling_suffix[index]+'.dat'
        #     UMI_dir=UMI_dir_base+sampling_suffix[index]+'/'
        #     os.system('python get_UMI_count_matrices.py -i '+ UMI_dir +' -t '+UMI_file+' -d '+UMI_distribution_file+' -g '+UMI_gene_file)
            
        # print('Getting pairwise distances between UMI matrices...')
        # UMI_distance_base='./Zeisel_UMI_gene_pairwise_SJ_subsample'
        # for index in range(1):
        #     UMI_distance_file=UMI_distance_base+sampling_suffix[index]+'.dat'
        #     UMI_gene_file=UMI_gene_file_base+sampling_suffix[index]+'.dat'
        #     os.system('python get_pairwise_distances.py '+UMI_gene_file+' '+UMI_distance_file+' '+str(num_proc))
            
        # print('Running eXpress....')
        # os.system('python run_express.py -r '+ref_transcriptome+' -n '+str(num_proc))

        # print('Processing eXpress cells....')
        # os.system('python process_xprs.py -n '+str(num_proc))
        # os.system('python t3i_to_expression_matrix.py')
            
        # print('Getting pairwise distances for eXpress...')
        # eXpress_gene_base='./Zeisel_express_gene_subsample'
        # eXpress_distance_base='./Zeisel_express_gene_distance_subsample'
        # for index in range(1):
        #     eXpress_distance_file=eXpress_distance_base+sampling_suffix[index]+'.dat'
        #     eXpress_gene_file=eXpress_gene_base+sampling_suffix[index]+'.dat'
        #     os.system('python get_pairwise_distances.py '+eXpress_gene_file+' '+eXpress_distance_file+' '+str(num_proc))
