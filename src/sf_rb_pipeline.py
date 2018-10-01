#!/usr/bin/python

# argv[1,2] = R1/R2.fastq files

import subprocess
import sys
import re
import os

steps_to_process = ["picard-tools", "GATK","Annovar"] # 
num_threads = "6"
# *********************************************************************
# DEFINITION OF PATHS

bwa_mem = ["bwa","mem"]
SamFormatConverter = ["java", "-jar", "/usr/share/java/picard.jar", "SamFormatConverter"]
SortSam  = ["java", "-jar", "/usr/share/java/picard.jar", "SortSam"]
MarkDuplicates = ["java", "-jar", "/usr/share/java/picard.jar", "MarkDuplicates"]
CollectHsMetrics = ["java", "-jar", "/usr/share/java/picard.jar", "CollectHsMetrics"]
BuildBamIndex = ["java", "-jar", "/usr/share/java/picard.jar", "BuildBamIndex"]
RealignerTargetCreator = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T RealignerTargetCreator"]
IndelRealigner = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T IndelRealigner"]
BaseRecalibrator = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T BaseRecalibrator"]
PrintReads = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T PrintReads"]
AnalyzeCovariates = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T AnalyzeCovariates"]
VariantFiltration = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T VariantFiltration"]
min_base_quality = "15"
min_read_length  = "25"
sliding_window   = "6:20" # window_size:min_average_quality
adapter_sequence = "AAGCAGTGGTATCAA"

fastq_r1_location = sys.argv[1]
fastq_r2_location = sys.argv[2]
reference_genome = sys.argv[3]


r1_base_filename = fastq_r1_location.replace(".fastq","").split("/")[-1]
r2_base_filename = fastq_r2_location.replace(".fastq","").split("/")[-1]
fastq_r1 = fastq_r1_location.split("/")[-1]
fastq_r2 = fastq_r2_location.split("/")[-1]
RBID = r1_base_filename.split("_")[0]
indels_vcf = "/media/thor/storage/rb_pipeline/mills", "/media/thor/storage/rb_pipeline/1000g"
output_directory = "/media/sf_RB_exome_project/"


# /DEFINITION OF PATHS
# *********************************************************************

# FASTQC
# =====================================================================
""" if("fastqc" in steps_to_process):
	print "running fastqc on the trimmed file"
	cmd = ["fastqc",fastq_r1_trimmed,fastq_r2_trimmed,fastq_r1_cutadapt,fastq_r2_cutadapt]
	print " ".join(cmd)
	subprocess.call(cmd)
	del cmd
"""

# BWA ALIGNER
# =====================================================================
#BWA MEM -aM -t 6 -R '@RG' reference_genome.fa fastq_r1 fastq_r2 > file.sam 
bwa_output_dir = output_directory +"bwa/"
bwa_sam= r1_base_filename + ".sam"
if("bwa-aligner" in steps_to_process):
	print "beginning alignment"
	print fastq_r1, fastq_r2
	out_file = bwa_output_dir +RBID+"_bwa.log"
	cmd = bwa_mem[:]
	cmd.append("-M")
	cmd.append("-t")
	cmd.append("6")
	#cmd.append(log)
	cmd.append(reference_genome)
	cmd.append(fastq_r1) 
	cmd.append(fastq_r2)
	cmd.append(">")
	cmd.append(bwa_sam)
		#print cmd
	print " ".join(cmd)
	p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
	print p
while True:
	out = p.stderr.read(1)
	if out == '' and p.poll() != None:
		break
	if out != '':
		sys.stdout.write(out)
		sys.stdout.flush()

# SAMFORMATCONVERTER
# =====================================================================
#SamFormatConverter I=bwa_output.sam O=converted_bwa_output.bam 
bwa_bam = r1_base_filename + ".bam"
if("picard-tools" in steps_to_process)
	cmd = SamFormatConverter[:]
	cmd.append("I=")
	cmd.append(bwa_sam)
	cmd.append("O=")
	cmd.append(bwa_output_dir,"/",bwa_bam)
	q = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
	print q

"""
# SORTSAM
# =====================================================================
#SortSam I=input.bam O=sorted.bam SORT_ORDER=coordinate 
picard_output_dir = output_directory +"picard/"
sortsam_bam= picard_output_dir+ r1_base_filename + ".sorted.bam"
if("picard-tools" in steps_to_process):
	subprocess.call(["mkdir",picard_output_dir])
	out_file = output_directory +"sortsam/"+RBID+"_sortsam.log"
	print "sortsam_bam:",sortsam_bam
	print "out_file: ",out_file
	cmd = SortSam[:]
	cmd.append("I="
	cmd.append(bwa_bam)
	cmd.append("O="
	cmd.apend(sortsam_bam)
	cmd.append("SORT_ORDER=coordinate"]
	print " ".join(cmd)
	
# MARKDUPLICATES
# =====================================================================
#MarkDuplicates I=input.bam O=marked_duplicates.bam M=marked_dup_metrics.txt 
marked_duplicates_bam= picard_output_dir+ r1_base_filename + ".marked_duplicates.bam"
marked_dup_metrics = r1_base_filename+"marked_dup_metrics.txt"
if("picard-tools" in steps_to_process):
	cmd = MarkDuplicates[:]
	cmd.append("I=")
	cmd.append(sortsam_bam)
	cmd.append("O="
	cmd.append(marked_duplicates_bam
	cmd.append("M=")
	cmd.append(marked_dup_metrics)
	print " ".join(cmd)
	with open(out_file,"w")as outf:
	 	subprocess.call(cmd,stdout=outf)

# COLLECTHSMETRICS
# =====================================================================
#CollectHsMetrics I=input.bam O=hsmetrics.txt R=reference_seq.fa BAIT_INTERVALS=bait_intervals.intervals TARGET_INTERVALS=target_intervals.intervals 
hs_metrics= picard_output_dir+ r1_base_filename + ".hsmetrics.txt"
if("picard-tools" in steps_to_process):
	cmd = [CollectHsMetrics]
	cmd += ["I=", marked_duplicates_bam]
	cmd += ["O=", r1_base_filename, ".hsmetrics.txt"]
	cmd += ["R=", reference_genome]
	cmd += ["BAIT_INTERVALS=", bait_intervals]
	cmd += ["TARGET_INTERVALS=", target_intervals]
	print " ".join(cmd)
	with open(out_file,"w")as outf:
	 	subprocess.call(cmd,stdout=outf)

# BUILDBAMINDEX
# =====================================================================
#BuildBamIndex I=input.bam 
if("picard-tools" in steps_to_process):
	cmd = [BuildBamIndex]
	cmd += ["I=", marked__duplicates_bam]
	print " ".join(cmd)
	with open(out_file,"w")as outf:
		subprocess.call(cmd,stdout=outf)

# REALIGNERTARGETCREATOR
# =====================================================================
#RealignerTargetCreator -T=RealignerTargetCreator -R=reference.fasta -I=input.bam --known=indels.vcf -o=forIndelRealigner.intervals
intervals_for_IR = RBID+"_"+"forIndelRealigner.intervals"
if("GATK" in steps_to_process):
	cmd = RealignerTargetCreator[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I")
	cmd.append(marked_duplicates_bam)
	cmd.append("-known")
	cmd.append(indels_vcf)
	cmd.append("-o")
	cmd.append(intervals_for_IR)
	
# INDELREALIGNER
# =====================================================================
#IndelRealigner -T=IndelRealigner -R=reference.fasta -I=input.bam -known=indels.vcf -targetIntervals=intervalListFromRTC.intervals -o=realignedBam.bam
realigned_bam = RBID+"realignedBam.bam"
if("GATK" in steps_to_process):
	cmd = IndelRealigner[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I")
	cmd.append(marked_duplicates_bam)
	cmd.append("-known")
	cmd.append(indels_vcf)
	cmd.append("-targetIntervals")
	cmd.append(intervals_for_IR)
	cmd.append("-o")
	cmd.append(realigned_bam)

# BASERECALIBRATOR
# =====================================================================
#BaseRecalibrator -T=IndelRealigner -R=reference.fasta -I=input.bam -known=indels.vcf -targetIntervals=intervalListFromRTC.intervals -o=realignedBam.bam 
recal_report = RBID+"recal_report.grp"
if("GATK" in steps_to_process):
	cmd = BaseRecalibrator[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I")
	cmd.append(realigned_bam)
	cmd.append("-knownSites")
	cmd.append(dbsnp_vcf)
	cmd.append("-o")
	cmd.append(recal_report)
	
# PRINTREADS
# =====================================================================
#PrintReads -T=PrintReads -R=reference.fasta -I=input.bam -BQSR= recalibration_report.grp -o=output.bam
recal_bam = RBID+"_recalibrated.bam"
if("GATK" in steps_to_process):
	cmd = PrintReads[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I")
	cmd.append(realigned_bam)
	cmd.append("-BQSR")
	cmd.append(recal_report)
	cmd.append("-o")
	cmd.append(recal_bam)
	
# ANALYZECOVARIATES
# =====================================================================
#AnalyzeCovariates -T=AnalyzeCovariates -R=reference.fasta -I=input.bam -BQSR= recalibration_report.grp -o=output.bam
recal_plot = RBID+"_recalQC.pdf"
if("GATK" in steps_to_process):
	cmd = AnalyzeCovariates[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-before")
	cmd.append(realigned_bam)
	cmd.append("-after")
	cmd.append(recal_report)
	cmd.append("-plots")
	cmd.append(recal_plot)
	
# VARIANTFILTRATION
# =====================================================================
#VariantFiltration -T=VariantFiltration -R=reference.fasta -o=output.vcf --variant= input.vcf --filterExpression= "AB < 0.2 || MQ0 > 50" --filterName= "SomeFilterName"
output_vcf = RBID+"_output.vcf"
input_vcf=
if("GATK" in steps_to_process):
	cmd = VariantFiltration[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-o")
	cmd.append(output_vcf)
	cmd.append("--variant")
	cmd.append(input_vcf)
	cmd.append("-o")
	cmd.append(recal_bam)
"""
