# to get the star annotation for each reads specified
grep 'NS500551:319:HKT73BGX2:1:13102:20662:10892\|NS500551:319:HKT73BGX2:1:13111:18829:16324\|NS500551:319:HKT73BGX2:1:21104:5233:16366\|NS500551:319:HKT73BGX2:1:22109:20474:7119\|NS500551:319:HKT73BGX2:2:11303:12836:6706\|NS500551:319:HKT73BGX2:3:13406:5445:5155\|NS500551:319:HKT73BGX2:4:22609:20034:19186\|NS500551:341:HHN7NBGX3:4:21409:13043:8697\|NS500551:341:HHN7NBGX3:4:23504:10329:7356\|NS500551:354:HHNKFBGX3:2:11311:10541:3101\|NS500551:354:HHNKFBGX3:3:11608:17710:15667' /home/chelseaju/DropSeq/Star/E31/exon_filter_read2gene.txt > ./AAATTCCTACTCAAGAGGGG_star.fastq

#find the read seq given read names as a file (don't > to a file!)
grep -f read_names/AAATTCCTACTCAAGAGGGG *reads_filtered.fastq -A 2

#grep all read names in a specified cluster (with size > 1)
grep "AAATTCCTACTCAAGAGGGG" All_cluster_other.fastq

makeblastdb -in ./blastdb/Mus_musculus_**_sequence.fa -dbtype nucl

blastn -task blastn-short -evalue 100 -query ./TTTGCGTCATAGCCGCGCTT.fa -db ../blast/*.fa

blastn -query TTTGCGTCATAGCCGCGCTT.fa -db ../blastdb/Mus_musculus_Large1_sequence.fa -task blastn-short -evalue 100

