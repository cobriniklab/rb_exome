#!/bin/bash

# argv[1] = directory with the input fastq files
# argv[2] = reference genome

list_to_process=$1 #"list of input fastq files"
reference_genome=$2
step_to_process=$3

for b in `cat $list_to_process`;do 
	j=`echo $b | sed 's/_1/_2/'`
	echo $i $j $reference_genome
	#fastqc $i $j
	./rb_batch_pipeline.py $b $j $reference_genome $step_to_process 
done
