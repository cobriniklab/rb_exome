#!/bin/bash

# argv[1] = directory with the input fastq files
# argv[2] = reference genome

list_to_process=$1 #"list of input files"

for b in `cat $list_to_process`;do 
	j=`echo $b | sed 's/.fastq.gz/_recalibrated.bam/'`
	echo $b $j
	fullpath = "/media/thor/storage/rb_pipeline/gatk/"+$j
	scp -r  $fullpath skevin@10.134.4.105:/dataVolume/storage/rb_pipeline/gatk/
done
