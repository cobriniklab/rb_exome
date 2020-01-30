#!/bin/bash -i

# argv[1] = bam file to process
# argv[2] = reference genome
# argv[3] = bin_size

# find output -type f \( -name "*recalibrated.bam" ! -name "*N_1_recalibrated.bam" \)

if [[ $# -eq 0 ]] ; then
    echo 'need to provide additional arguments:
        argv[1] = bam file
        argv[2] = reference genome
        argv[3] = bin_size
        argv[4] = matched normals present'
        exit 0
fi

bam_in=$1	#list of input bam files"
reference_genome=$2
bin_size=$3
matched_normals=$4

# configure files
# get path and prefix for input bam file
sample_name=`basename $bam_in | sed 's;_.*;;'`
sample_extracted=`echo "$sample_name"_extracted.txt`
# define output for sample and reference mpileups
sample_mpileup=`echo "output/bafsegmentation/data/"$sample_name".mpileup"`
# define output for mpileup2baf for sample and reference
sample_txt=`echo "output/bafsegmentation/extracted/"$sample_name".txt"`
sample_igv=`echo "output/bafsegmentation/data/"$sample_name"_read_counts.igv"`

# configure normals
# get path and prefix for matched-normal bam file
reference_path=`echo $bam_in | sed 's/T/N/g; s/CL/N/g'`
echo $reference_path
reference_name=`echo $reference_path | sed 's/.*\///' | sed 's;_.*;;'`
reference_extracted=`echo "$reference_name"_extracted.txt`
reference_mpileup=`echo "output/bafsegmentation/data/"$reference_name".mpileup"`
reference_txt=`echo "output/bafsegmentation/extracted/"$reference_name".txt"`

# tumor/normal comparison
normal_sample_names_txt="/home/skevin/rb_pipeline/output/bafsegmentation/extracted/normal_sample_names.txt"
echo -e "FilenameAssay"'\t'"FilenameNormal" > $normal_sample_names_txt

echo $bam_in

echo -e "$sample_extracted"'\t'"$reference_extracted" >> $normal_sample_names_txt

# sample_igv=`echo "output/bafsegmentation/data/log2-"$sample_name"_read_counts.igv"`
sample_igv=`echo "output/bafsegmentation/data/"$sample_name"_read_counts.igv"`

echo $bam_in $sample_name $sample_mpileup $sample_txt
echo $reference_path $reference_name $reference_mpileup $reference_txt 

# run samtools on each sample and matching reference
if [ ! -f $sample_mpileup ]; then
	samtools mpileup --min-MQ 35 --min-BQ 20 -f $reference_genome -l output/bafsegmentation/data/all_samples_merged2.bed $bam_in | awk '$4 > 20' > $sample_mpileup 
fi
if [ ! -f $reference_mpileup ]; then
	samtools mpileup --min-MQ 35 --min-BQ 20 -f $reference_genome -l output/bafsegmentation/data/all_samples_merged2.bed $reference_path | awk '$4 > 20' > $reference_mpileup 
fi

#~ # run mpileup2baf on each sample
perl bin/mpileup2baf.pl < $sample_mpileup > $sample_txt
perl bin/mpileup2baf.pl < $reference_mpileup > $reference_txt

echo "add logr ratios from "$sample_igv" to "$sample_txt" before running bafsegmentation"
echo ~/rb_pipeline/src/LOH/add_log_R_value_to_BAF.py $sample_igv $sample_txt $reference_txt $bin_size
~/rb_pipeline/src/LOH/add_log_R_value_to_BAF.py $sample_igv $sample_txt $reference_txt $bin_size

#~ split samples for bafsegmentation
r_log_file=`echo $sample_txt | sed 's/.txt/_with_Rlog.csv/'`
echo 
perl /home/skevin/rb_pipeline/output/bafsegmentation/split_samples.pl --data_file=$r_log_file --output_directory=output/bafsegmentation/extracted

#~ run bafsegmentation
perl output/bafsegmentation/BAF_segment_samples.pl --matched_normals="$matched_normals"
mkdir output/bafsegmentation/segmented/$sample_name/
cp output/bafsegmentation/segmented/AI_regions.txt output/bafsegmentation/segmented/$sample_name/AI_regions.txt

