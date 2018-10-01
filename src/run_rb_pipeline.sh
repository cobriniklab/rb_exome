#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo 'need to provide additional arguments:
        argv[1] = directory with the input fastq files
        argv[2] = organism
        argv[3] = output dir
        argv[4] = step to process'
    exit 0
fi

directory_with_fastq_files=$1 #"input_fastq"
organism=$2
outdir=$3
step_to_process=$4

for i in `find $directory_with_fastq_files -name "*_1.fastq.gz" | sort -V`;do # !!!! REMOVE tail -n +2 !!!
	j=`echo $i | sed 's/_1/_2/'`
	echo $i $j $reference_genome
	#fastqc $i $j
	./src/rb_pipeline.py -1 $i -2 $j -d $outdir -rg "/dataVolume/storage/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa" -s $step_to_process
done
