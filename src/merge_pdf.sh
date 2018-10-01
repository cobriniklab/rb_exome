#!/bin/bash
#24-T_vs_N-all_chrom.pdf 

input_directory=$1

for T in `find $input_directory -wholename "*T*N*/all_chrom.pdf"`; do
        C=`echo $T | sed 's/T/CL/'`
        file=`echo $T | sed 's/.*\///'`
        T_name=`echo $T | cut -c 8-11 -n`
	outpdf=""$T_name"_and_CL_merged_all_chrom.pdf"
        echo $outpdf
	pdftk $T $C cat output $outpdf
        pdfjam $T $C --nup 1x2 --landscape --outfile $outpdf
done

pdftk *_merged_all_chrom.pdf cat output "all_chrom_total.pdf"

~                                                     
