#!/bin/bash -i

# argv[1] = bam file to process
# argv[2] = reference genome
# argv[3] = bin_size

if [[ $# -eq 0 ]] ; then
    echo 'need to provide additional arguments:
        argv[1] = list of bam files in experiment
        argv[2] = reference genome
        argv[3] = bin_size'
    exit 0
fi

list_bam_in=$1	#list of input bam files"
reference_genome=$2
bin_size=$3

normal_sample_names_txt="/home/skevin/rb_pipeline/output/bafsegmentation/extracted/normal_sample_names.txt"
echo -e "FilenameAssay"'\t'"FilenameNormal" > $normal_sample_names_txt 

for i in `cat $list_bam_in`;do
	echo $i
	# get path and prefix for input bam file
	sample_name=`echo $i | sed 's/.*\///' | sed 's;_.*;;'`

	# get path and prefix for matched-normal bam file
	reference_path=`echo $i | sed 's/T/N/g; s/CL/N/g'`
	echo $reference_path
	reference_name=`echo $reference_path | sed 's/.*\///' | sed 's;_.*;;'`

	# define output for sample and reference mpileups
	sample_mpileup=`echo "output/bafsegmentation/data/"$sample_name".mpileup"`
	reference_mpileup=`echo "output/bafsegmentation/data/"$reference_name".mpileup"`

	# define output for mpileup2baf for sample and reference
	sample_txt=`echo "output/bafsegmentation/extracted/"$sample_name".txt"`
	reference_txt=`echo "output/bafsegmentation/extracted/"$reference_name".txt"`
	
	sample_igv=`echo "output/bafsegmentation/data/log2-"$sample_name"_read_counts.igv"`

	echo $i $sample_name $sample_mpileup $sample_txt
	echo $reference_path $reference_name $reference_mpileup $reference_txt 
	
	# run samtools on each sample and matching reference
	if [ ! -f $sample_mpileup ]; then
		samtools mpileup --min-MQ 35 --min-BQ 20 -f $reference_genome -l output/bafsegmentation/data/all_samples_merged.bed $i | awk '$4 > 20' > $sample_mpileup 
	fi
	if [ ! -f $reference_mpileup ]; then
		samtools mpileup --min-MQ 35 --min-BQ 20 -f $reference_genome -l output/bafsegmentation/data/all_samples_merged.bed $reference_path | awk '$4 > 20' > $reference_mpileup 
	fi

	#~ # run mpileup2baf on each sample
	perl bin/mpileup2baf.pl < $sample_mpileup > $sample_txt
	perl bin/mpileup2baf.pl < $reference_mpileup > $reference_txt
	
	#~ # add logr ratios to baf.txt file before running bafsegmentation
	~/rb_pipeline/src/add_log_R_value_to_BAF.py $sample_igv $sample_txt $reference_txt $bin_size
	
	#~ split samples for bafsegmentation
	r_log_file=`echo $sample_txt | sed 's/.txt/_with_Rlog.csv/'`
	perl /home/skevin/rb_pipeline/output/bafsegmentation/split_samples.pl --data_file=$r_log_file --output_directory=output/bafsegmentation/extracted
	
	#~ run bafsegmentation
	perl output/bafsegmentation/BAF_segment_samples.pl
	mkdir output/bafsegmentation/segmented/$sample_name/
	cp output/bafsegmentation/segmented/AI_regions.txt output/bafsegmentation/segmented/$sample_name/AI_regions.txt
done 
	#~ # run bafsegmentation
	#~ baf_input=`echo $sample_txt | sed 's/.txt/_with_Rlog.csv/'`
	#~ echo $baf_input
	#~ mkdir output/bafsegmentation/segmented/$sample_name
	#~ perl output/bafsegmentation/split_samples.pl --data_file=$baf_input --output_directory=output/bafsegmentation/extracted
	#~ perl output/bafsegmentation/BAF_segment_samples.pl --input_directory=output/bafsegmentation/extracted \
	#~ --output_directory=output/bafsegmentation/segmented/$sample_name/ \
	#~ --plot_directory=output/bafsegmentation/plots --ai_size=5 \
	#~ --matched_normals=TRUE --run_matched_normals=FALSE
	
	#~ #concatentat AI regions files
	#~ segment_out_dir="output/bafsegmentation/segmented"
	#~ echo $segment_out_dir
	#~ head -1 "$segment_out_dir"/"$sample_name"/*.txt > "$segment_out_dir"/all_AI_regions.txt; tail -n +2 -q "$segment_out_dir"/*/AI_regions.txt >> "$segment_out_dir"/all_AI_regions.txt
	
		
#~ done
