#!/usr/bin/bash

for line in $(cat vcf_list.txt);do # !!!! REMOVE tail -n +$
        j=`echo $i | sed 's/_1/_2/'`
        echo $i $j $reference_genome
        #fastqc $i $j
        ./rb_pipeline.py $i $j $reference_genome $step_to_process 
done


