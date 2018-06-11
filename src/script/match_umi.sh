#!/bin/sh
IFS=$'\r\n' GLOBIGNORE='*' command eval  'readnames=($(cat /home/zgy_ucla_cs/Research/DropSeq/data/read_names/E31_readnames.txt))'
#echo "${readnames[@]}"
echo "Done with loading readmanes"
for name in ${readnames[@]}
do
#echo $name
grep -m1 -A3 $name /home/zgy_ucla_cs/Research/DropSeq/data/Reps/E31_1.fastq
done
