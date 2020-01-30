#!/bin/bash

for i in `find ~/rb_pipeline/output/gatk/ -name "*N_1_recalibrated.bam"`;do i
    echo $i; 
    j=`basename $i | sed 's/_recalibrated.bam/all_readcount.txt/g'`; 
    bam-readcount -l ~/rb_pipeline/results/SNV/reported_m2_reynolds_somatic_variants.bed -f ~/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa $i > $j;
done

