'''
Sample run: 
python extract_barcode_bam.py ../../data/7_uniqueAligned_clean_gene_exon_tagged_sorted.bam ../../data/cell_names_filtered.txt /home/zgy_ucla_cs/Research/DropSeq/data/cells_temp

Reads:
@NS500551:319:HKT73BGX2:1:11101:14704:1049 2:N:0:TAAGGCGA
TTGTTTGATTGTTTGTTTTGTTTTGTTTTCTTTTGAAATAAAACAGGAAAG
+
A/AA/E/EAEAE<EA//<E6/E/E/E/E</E/<EA/6E/<E/E/EAEEE/A
'''
import regex as re
import sys
from argparse import ArgumentParser
from pysam import AlignmentFile
def extract_barcode(sam, barcode_file, outdir):

    # Create the hash set for cell names
    fin = open(barcode_file, 'r')
    barcodes_filtered = set()
    for line in fin:
        line = line.strip()
        barcodes_filtered.add(line)

    print(len(barcodes_filtered))

    sam_file = AlignmentFile(sam, mode='r')
    #filter_file = AlignmentFile("-", mode='wh', template=sam_file)
    track = sam_file.fetch(until_eof=True)
    for i, aln in enumerate(track):
        # if aln.is_unmapped:
            # continue
        # print(i)

        ''' Error to use query_alignment_sequence, use query_sequence instead? '''
        reads_name, reads, cell_barcode, umi, quality = aln.qname, aln.query_sequence, aln.get_tag('XC'), aln.get_tag('XM'), aln.qual
        # print(reads_name, reads, cell_barcode, umi, quality)
        # print(reads)
        # print(quality)
        
        if cell_barcode in barcodes_filtered:
            # print(reads_name, reads, cell_barcode, umi, quality)
            if len(reads) != len(aln.qual):
                print("Error, skipped:", reads, quality)
                continue
                
            fout_umi = open(outdir + '/' + cell_barcode + '.umi', 'a+')
            fout_umi.write(umi + '\n')

            fout_fq = open(outdir + '/' + cell_barcode + '.fastq', 'a+')
            fout_fq.write('@' + reads_name + '\n')
            fout_fq.write(reads + '\n')
            fout_fq.write('+\n')
            fout_fq.write(quality + '\n')
        if i % 100000 == 0:
            print(i/209400000.0)

if __name__ == "__main__":
    parser = ArgumentParser("extract reads/alignments from a single cell")
    parser.add_argument("file", help="A SAM or FASTQ file")
    parser.add_argument("barcodes", help="barcode of the cell to extract")
    parser.add_argument("outdir", help="output of the individual cells")
    args = parser.parse_args()
    extract_barcode(args.file, args.barcodes, args.outdir)
