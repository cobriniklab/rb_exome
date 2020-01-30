#!/usr/bin/python

import sys
import os.path

# argv[] = fastq R1 files to be checked
# FASTQ/14-CL_1.fastq.gz
summary_file = "/dataVolume/storage/rb_pipeline/srt_rb_summary.csv"
missings_file = "/dataVolume/storage/rb_pipeline/srt_rb_missing.csv"
header = [""]
header += ["R1.fastq","R2.fastq"]
header += ["R1_fastqc.html","R2_fastqc.html"]
header += [".bam","_sorted.bam"]
header += ["_marked_duplicates.bam"]
header += ["_hsmetrics.txt"]
header += ["_ismetrics.txt"]
header += ["_alignsummetrics.txt"]
header += ["_marked_duplicates.bai"]
header += ["base_recalibrator"]
header += ["recalibrated.bam"]
header += ["mutect2.vcf"]


with open(summary_file, 'w') as f, open(missings_file, 'w') as file:
	f.write("\t".join(header))	
	for f1 in sys.argv[1:]:
		f1 = f1.split("/")[-1]
		f2 = f1.replace("_1","_2")
		#print f1,f2
		expected_files = [f1,f2] # fastq R1
		b1 = f1.split(".")[0] # 14-CL_1
		b2 = f2.split(".")[0] # 14-CL_2
		expected_files.append("fastqc/"+b1+"_fastqc.html")
		expected_files.append("fastqc/"+b2+"_fastqc.html")
		expected_files.append("picard/"+b1+".bam")# picard/14-CL_1
		expected_files.append("picard/"+b1+"_sorted.bam")
		expected_files.append("picard/markduplicates/"+b1+"_marked_duplicates.bam")
		expected_files.append("picard/alignment_qc/collecthsmetrics/"+b1 + "_hsmetrics.txt")
		expected_files.append("picard/alignment_qc/collectinsertsizemetrics/"+b1 + "_ismetrics.txt")
		expected_files.append("picard/alignment_qc/collectalignmentsummarymetrics/"+b1 + "_alignsummetrics.txt")
		expected_files.append("picard/markduplicates/"+b1+ "_marked_duplicates.bam")
		expected_files.append("gatk/"+b1+ "_recal_report.table")
		expected_files.append("gatk/"+b1+ "_recalibrated.bam")
		expected_files.append("gatk/mutect2/"+b1+ "_mutect2.vcf")
	
		out = [f1]
		for i in expected_files:
			if(os.path.isfile(i)):
				out.append(str(os.path.getsize(i)))
			else:
				out.append("N")
		f.write("\t".join(out)+"\n")

		for i in expected_files:
                	if not(os.path.isfile(i)):
				file.write("\n"+os.path.basename(i))

