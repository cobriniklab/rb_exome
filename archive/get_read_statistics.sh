#!/bin/bash

# $1 = folder where to search for align_summary.txt
# $2 = output filename with alignment statistics

echo "cell	mapped	uniquly_mapped	concordant" > $2
for i in `find $1 -name "align_summary.txt" | sort -t '_' -k2g`;do
	summary_lines=`wc -l $i | cut -f 1 -d' '`
	name=`echo $i | sed 's/.*\(Hu_[^_]*\).*/\1/'`
	if [ $summary_lines -eq 14 ];then
		left_mapped=`sed -n '3p' $i | sed 's/.*Mapped *: *\([0-9]*\).*/\1/'`
		right_mapped=`sed -n '7p' $i | sed 's/.*Mapped *: *\([0-9]*\).*/\1/'`
		
		left_left_multiple_alignments=`sed -n '4p' $i | sed 's/.*of these *: *\([0-9]*\).*/\1/'`
		right_left_multiple_alignments=`sed -n '8p' $i | sed 's/.*of these *: *\([0-9]*\).*/\1/'`
		
		mapped=`echo $left_mapped + $right_mapped | bc`
		uniqly_mapped=`echo $mapped - $left_left_multiple_alignments - $right_left_multiple_alignments | bc`
		
		pairs_mapped=`sed -n '11p' $i | sed 's/.*Aligned pairs *: *\([0-9]*\).*/\1/'`
		pairs_multi_mapped=`sed -n '12p' $i | sed 's/.*of these *: *\([0-9]*\).*/\1/'`
		pairs_discordant=`sed -n '13p' $i | sed 's/ *\([0-9]*\).*/\1/'`
		
		concordant=`echo "( $pairs_mapped - $pairs_discordant) *2" | bc`
		
		# there are 4 lines per read, however R1 accounts only for half reads. x/4 *2 = x/2
		echo $name"	"$mapped"	"$uniqly_mapped"	"$concordant
	else
		(>&2 echo "error with" $i)
		echo $name"	0	0	0"
	fi
done >> $2
