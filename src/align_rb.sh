#!/bin/bash

# argv[1] = directory with the input fastq files
# argv[2] = reference genome

directory_with_fastq_files=$1 #"input_fastq"
reference_genome=$2

for i in `find $directory_with_fastq_files -name "*_1.fastq.gz"`;do # !!!! REMOVE tail -n +2 !!!
	j=`echo $i | sed 's/_1/_2/'`
	echo $i $j $reference_genome
	#fastqc $i $j
	r1_sam=`echo $i | sed 's/.fastq.gz/.sam/'`
	bwa mem -M -t 6 ucsc.hg19.fasta $i $j > $r1_sam
done
