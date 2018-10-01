#!/bin/bash -i

# argv[1] = bam file to process
# argv[2] = reference genome

if [[ $# -eq 0 ]] ; then
    echo 'need to provide additional arguments:
        argv[1] = list of all vcf files in experiment
        argv[2] = reference genome
		argv[3] = outdir' #output/bafsegmentation/data
    exit 0
fi

list_vcf_in=$1 #list of input vcf files"
reference_genome=$2
outdir=$3

for i in `cat $list_vcf_in`;do
	filename=`echo $i | sed 's/.*\///' | sed 's;_.*;;'`
	snv_bed=`echo $outdir"/"$filename"_snv.bed"`
	del_bed=`echo $outdir"/"$filename"_del.bed"`
	sample_bed=`echo $outdir"/"$filename"_final.bed"`
	echo $vcf_to_process $filename $snv_bed $del_bed 

	# create merged snp bed file
	vcf2bed --deletions < $i > $del_bed
	vcf2bed --snvs < $i > $sample_bed #snv_bed
	bedops --everything $del_bed $snv_bed > $sample_bed
done

bedops --everything $outdir/*_final.bed > $outdir/all_samples_merged.bed
