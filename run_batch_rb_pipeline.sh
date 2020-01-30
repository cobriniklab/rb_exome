#!/bin/bash
# argv[1] = directory with the input fastq files
# argv[2] = reference genome

if [[ $# -eq 0 ]] ; then
    echo 'need to provide additional arguments:
        argv[1] = directory with the input fastq files
        argv[2] = organism
        argv[3] = output dir
        argv[4] = step to process'
    exit 0
fi

list_to_process=$1 #"list of input fastq files"
organism=$2
outdir=$3
step_to_process=$4

for b in `cat $list_to_process | sort -V`;do
        j=`echo $b | sed 's/_1/_2/'`
        cell_name=`echo $b | sed 's#.*/\(.*\)_1.*#\1#'`
        echo $b $j $cell_name
        ./src/rb_pipeline.py -1 $b -2 $j -d $outdir/ -rg "/dataVolume/storage/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa" -s $step_to_process
done

if [[ $step_to_process == "create_pon" ]];
then
    # echo "sucess"
    ./src/SNV/create_pon.py $list_to_process $outdir
fi