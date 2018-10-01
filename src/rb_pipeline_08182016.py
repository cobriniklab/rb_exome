#!/usr/bin/python

# argv[1,2] = R1/R2.fastq files

import subprocess
from subprocess import PIPE
import sys
import re
import os

steps_to_process = ["bwa-aligner", "picard-tools", "GATK","Annovar"] # 
num_threads = "6"
# *********************************************************************
# DEFINITION OF PATHS

BwaMem = ["bwa", "mem"]
SamFormatConverter = ["java",  "-jar", "/usr/share/java/picard.jar", "SamFormatConverter"]
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

output_directory = "/media/thor/storage1/rb_pipeline/"
bwa_output_directory = output_directory
picard_output_dir = output_directory 
picard_tmp_dir = output_directory+"tmp/"
gatk_output_dir = output_directory + "gatk/"
r1_base_filename = fastq_r1_location.replace(".fastq.gz","").split("/")[-1]
r2_base_filename = fastq_r2_location.replace(".fastq.gz","").split("/")[-1]
fastq_r1 = fastq_r1_location.split("/")[-1]
fastq_r2 = fastq_r2_location.split("/")[-1]
RBID = r1_base_filename.split("_")[0]
bwa_sam= r1_base_filename + ".sam"

# /DEFINITION OF PATHS
# *********************************************************************

"""
# FASTQC
# =====================================================================
if("fastqc" in steps_to_process):
	print "running fastqc on the trimmed file"
	cmd = ["fastqc",fastq_r1_trimmed,fastq_r2_trimmed,fastq_r1_cutadapt,fastq_r2_cutadapt]
	print " ".join(cmd)
	subprocess.call(cmd)
	del cmd
"""

# BWA ALIGNER
# =====================================================================
#BWA MEM -aM -t 6 -R '@RG' reference_genome.fa fastq_r1 fastq_r2 > file.sam 
bwa_sam= r1_base_filename + ".sam"
if("bwa-aligner" in steps_to_process):
	print "beginning alignment"
	print fastq_r1, fastq_r2
	cmd = BwaMem[:]
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
	subprocess.call(cmd, shell=True) 
	del cmd
"""
# SAMFORMATCONVERTER
# =====================================================================
#SamFormatConverter I=bwa_output.sam O=converted_bwa_output.bamt
picard_bam = picard_output_dir+r1_base_filename + ".bam"
if("picard-tools" in steps_to_process):
	cmd = SamFormatConverter[:]
	cmd.append("I=")
	cmd.append(bwa_sam)
	cmd.append("O=")
	cmd.append(picard_bam)
	cmd.append("MAX_RECORDS_IN_RAM=5000000")
	cmd.append("TMP_DIR=")
	cmd.append(picard_tmp_dir)
	print " ".join(cmd)
	#sam_convert = subprocess.Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
	#sam_convert.communicate()
	subprocess.call(cmd)
	#with open (bwa_sam, "r") as i:
		#sam_convert = subprocess.Popen(cmd, shell=True, stdin=i, stdout=PIPE, stderr=PIPE)
		#sam_convert.stdout = open(picard_bam, 'w')
		#sam_convert.wait()
		#outbam = sam_convert.stdout.read()
		#with open(picard_bam, "w") as o:
		#	o.write(outbam)
	del cmd

# SORTSAM
# =====================================================================
#SortSam I=input.bam O=sorted.bam SORT_ORDER=coordinate 
sortsam_bam = picard_output_dir+ r1_base_filename + "_sorted.bam"
if("picard-tools" in steps_to_process):
	subprocess.call(["mkdir",picard_output_dir])
	print "sortsam_bam:",sortsam_bam
	cmd = SortSam[:]
	cmd.append("I=")
	cmd.append(picard_bam)
	cmd.append("O=")
	cmd.append(sortsam_bam)
	cmd.append("SORT_ORDER=coordinate")
	print " ".join(cmd)
	sort_sam = subprocess.Popen(cmd, shell=True, stdin=sam_convert.stdout, stdout=subprocess.PIPE, stderr=None)
	while True:
		out = sort_sam.stderr.read(1)
		if out == '' and sort_sam.poll() != None:
			break
		if out != '':
			sys.stdout.write(out)
			sys.stdout.flush()
	del cmd

# MARKDUPLICATES
# =====================================================================
#MarkDuplicates I=input.bam O=marked_duplicates.bam M=marked_dup_metrics.txt
marked_duplicates_bam= picard_output_dir+ r1_base_filename + "_marked_duplicates.bam"
marked_dup_metrics = picard_output_dir+r1_base_filename+"_marked_dup_metrics.txt"
picard_tmp_dir = picard_output_dir+"tmp/"
if("picard-tools" in steps_to_process):
	cmd = MarkDuplicates[:]
	cmd.append("I=")
	cmd.append(sortsam_bam)
	cmd.append("O=")
	cmd.append(marked_duplicates_bam)
	cmd.append("M=")
	cmd.append(marked_dup_metrics)
	cmd.append("TMP_DIR=")
	cmd.append(picard_tmp_dir)
	print " ".join(cmd)
	mark_dup = subprocess.Popen(cmd, shell=True, stdin=sort_sam.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	while True:
		out = mark_dup.stderr.read(1)
		if out == '' and mark_dup.poll() != None:
			break
		if out != '':
			sys.stdout.write(out)
			sys.stdout.flush()
	del cmd

# COLLECTHSMETRICS
# =====================================================================
#CollectHsMetrics I=input.bam O=hsmetrics.txt R=reference_seq.fa BAIT_INTERVALS=bait_intervals.intervals TARGET_INTERVALS=target_intervals.intervals 
hs_metrics= output_directory+ r1_base_filename + "_hsmetrics.txt"
bait_intervals = output_directory+"agilent_baits_hg19.interval_list"
target_intervals = output_directory+"agilent_targets_hg19.interval_list"
if("picard-tools" in steps_to_process):
	cmd = CollectHsMetrics[:]
	cmd.append("I=")
	cmd.append(marked_duplicates_bam)
	cmd.append("O=") 
	cmd.append(hs_metrics)
	cmd.append("R=")
	cmd.append(reference_genome)
	cmd.append("BAIT_INTERVALS=")
	cmd.append(bait_intervals)
	cmd.append("TARGET_INTERVALS=")
	cmd.append(target_intervals)
	print " ".join(cmd)
	collect_hs_met = subprocess.Popen(cmd, shell=True, stdin=mark_dup.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	while True:
		out = collect_hs_met.stderr.read(1)
		if out == '' and collect_hs_met.poll() != None:
			break
		if out != '':
			sys.stdout.write(out)
			sys.stdout.flush()
	del cmd



# BUILDBAMINDEX
# =====================================================================
#BuildBamIndex I=input.bam 
if("picard-tools" in steps_to_process):
	cmd = BuildBamIndex[:]
	cmd.append("I=")
	cmd.append(marked_duplicates_bam)
	print " ".join(cmd)
	build_bai = subprocess.Popen(cmd, shell=True, stdin=collect_hs_met.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	while True:
		out = build_bai.stderr.read(1)
		if out == '' and build_bai.poll() != None:
			break
		if out != '':
			sys.stdout.write(out)
			sys.stdout.flush()
	del cmd

# REALIGNERTARGETCREATOR
# =====================================================================
#RealignerTargetCreator -T=RealignerTargetCreator -R=reference.fasta -I=input.bam --known=indels.vcf -o=forIndelRealigner.intervals
intervals_for_IR = gatk_output_dir+r1_base_filename+"_forIndelRealigner.intervals"
indels_vcf_1 = output_directory+"Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
indels_vcf_2 = output_directory+"1000G_phase1.indels.hg19.sites.vcf"
if("GATK" in steps_to_process):
	cmd = RealignerTargetCreator[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I")
	cmd.append(marked_duplicates_bam)
	cmd.append("-known")
	cmd.append(indels_vcf_1)
	cmd.append("-known")
	cmd.append(indels_vcf_2)
	cmd.append("-o")
	cmd.append(intervals_for_IR)
	print " ".join(cmd)
	realign_tc = subprocess.Popen(cmd, shell=True, stdin=build_bai.stdout, stderr=subprocess.PIPE)
	while True:
		out = realign_tc.stderr.read(1)
		if out == '' and realign_tc.poll() != None:
			break
		if out != '':
			sys.stdout.write(out)
			sys.stdout.flush()
	del cmd
	
# INDELREALIGNER
# =====================================================================
#IndelRealigner -T=IndelRealigner -R=reference.fasta -I=input.bam -known=indels.vcf -targetIntervals=intervalListFromRTC.intervals -o=realignedBam.bam
realigned_bam = gatk_output_dir+r1_base_filename+"_realigned.bam"
if("GATK" in steps_to_process):
	cmd = IndelRealigner[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I")
	cmd.append(marked_duplicates_bam)
	cmd.append("-known")
	cmd.append(indels_vcf_1)
	cmd.append("-known")
	cmd.append(indels_vcf_2)
	cmd.append("-targetIntervals")
	cmd.append(intervals_for_IR)
	cmd.append("-o")
	cmd.append(realigned_bam)
	print " ".join(cmd)
	indel_realign = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
	while True:
		out = indel_realign.stderr.read(1)
		if out == '' and indel_realign.poll() != None:
			break
		if out != '':
			sys.stdout.write(out)
			sys.stdout.flush()
	del cmd

# BASERECALIBRATOR
# =====================================================================
#BaseRecalibrator -T=IndelRealigner -R=reference.fasta -I=input.bam -known=indels.vcf -targetIntervals=intervalListFromRTC.intervals -o=realignedBam.bam 
recal_report = gatk_output_dir+r1_base_filename+"_recal_report.grp"
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
recal_bam = gatk_output_dir+r1_base_filename+"_recalibrated.bam"
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


