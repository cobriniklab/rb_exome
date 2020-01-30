#!/usr/bin/python

# argv[1,2] = R1/R2.fastq files
# argv[3] = reference genome (*.fasta)
# argv[4] = step to process

import sys
import re
import os
import subprocess
import gzip


num_threads = "6"
# *********************************************************************
# DEFINITION OF PATHS

BwaMem = ["bwa", "mem"]
SamFormatConverter = ["java",  "-jar", "/usr/share/java/picard.jar", "SamFormatConverter"]
SortSam  = ["java", "-jar", "/usr/share/java/picard.jar", "SortSam"]
MarkDuplicates = ["java", "-jar", "/usr/share/java/picard.jar", "MarkDuplicates"]
CollectHsMetrics = ["java", "-jar", "/usr/share/java/picard.jar", "CollectHsMetrics"]
CollectIsMetrics = ["java", "-jar", "/usr/share/java/picard.jar", "CollectInsertSizeMetrics"]
CollectAlignmentSummaryMetrics = ["java", "-jar", "/usr/share/java/picard.jar", "CollectAlignmentSummaryMetrics"]
BuildBamIndex = ["java", "-jar", "/usr/share/java/picard.jar", "BuildBamIndex"]
RealignerTargetCreator = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "RealignerTargetCreator"]
IndelRealigner = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "IndelRealigner"]
BaseRecalibrator = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "BaseRecalibrator"]
PrintReads = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "PrintReads"]
DepthOfCoverage = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "DepthOfCoverage"]
AnalyzeCovariates = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "AnalyzeCovariates"]
VariantFiltration = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "VariantFiltration"]
SelectVariants = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "SelectVariants"]
Mutect2 = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "MuTect2"]
Mutect = ["java", "-jar", "/usr/share/java/muTect-1.1.4.jar"]
HaplotypeCaller = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "HaplotypeCaller"]
GenotypeGVCFs = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "GenotypeGVCFs"]
CopyWriter = ["/dataVolume/storage/rb_pipeline/copywriter/copywriter_cobrinik.R"]
HiResCopyWriter = ["/dataVolume/storage/rb_pipeline/hi_res_copywriter/hi_res_copywriter_cobrinik.R"]


min_base_quality = "15"
min_read_length  = "25"
sliding_window   = "6:20" # window_size:min_average_quality
adapter_sequence = "AAGCAGTGGTATCAA"

fastq_r1_location = sys.argv[1]
fastq_r2_location = sys.argv[2]
reference_genome = sys.argv[3]
steps_to_process = sys.argv[4].split(",")

output_directory = "/dataVolume/storage/rb_pipeline/"
storage_directory = "/dataVolume/storage/rb_pipeline/"
fastqc_output_dir = storage_directory+"fastqc/"
fastqc_storage_dir = storage_directory+"fastqc/"
bwa_storage_dir = storage_directory+"bwa/"
bwa_storage_dir = storage_directory+"bwa/"
picard_storage_dir = storage_directory+"picard/"
picard_storage_dir = storage_directory+"picard/"
picard_tmp_dir = storage_directory+"tmp/"
alignment_qc_storage_dir = storage_directory+"picard/alignment_qc/"
gatk_storage_dir = storage_directory + "gatk/"
gatk_storage_dir = storage_directory + "gatk/"
mutect2_storage_dir=storage_directory+"gatk/mutect2/"
mutect2_storage_dir=storage_directory+"gatk/mutect2/"
mutect_storage_dir=storage_directory+"gatk/mutect/"
haplocaller_storage_dir=storage_directory+"gatk/haplocaller/"
annovar_storage_dir = storage_directory + "annovar/"
r1_base_filename = fastq_r1_location.replace(".fastq.gz","").split("/")[-1]
r2_base_filename = fastq_r2_location.replace(".fastq.gz","").split("/")[-1]
fastq_r1 = fastq_r1_location.split("/")[-1]
fastq_r2 = fastq_r2_location.split("/")[-1]
RBID = r1_base_filename.split("_")[0]
bwa_sam= r1_base_filename + ".sam"
picard_bam = storage_directory+r1_base_filename + ".bam"
sortsam_bam = picard_storage_dir+ r1_base_filename + "_sorted.bam"

# /DEFINITION OF PATHS
# *********************************************************************

# READ GROUP PARSING
# =====================================================================
if("rg_parse" in steps_to_process) or ("all" in steps_to_process):
	with gzip.open(fastq_r1_location) as fasta:
		s = fasta.readline().split(":")
		RGID =  "_".join(s[2:5])
		RGSM = r1_base_filename
		RGPL = "Illumina"
		read_groups = "@RG\tID:"+RGID+"\tSM:"+RGSM+"\tPL:"+RGPL

# FASTQC
# =====================================================================
fastqc_storage_dir = storage_directory + "fastqc/"
fastqc_storage = fastqc_storage_dir+r1_base_filename+"_fasqc.zip"
if not os.path.isdir(fastqc_storage_dir):
	subprocess.call(["mkdir", fastqc_storage_dir])
if(not os.path.isfile(fastqc_storage)) and(("fastqc" in steps_to_process) or ("all" in steps_to_process)):
	cmd = ["fastqc",fastq_r1_location,fastq_r2_location]
	print " ".join(cmd)
	out_file = fastqc_storage_dir + "fastqc.log"
	err_file = fastqc_storage_dir + "fastqc.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
		 	errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
			if(errcode == 0):
				print "fastqc finished successfully"
			else:
				print "fastqc failed !!!!"
				del steps_to_process[:]
	del cmd

# BWA ALIGNER
# =====================================================================
#BWA MEM -aM -t 6 -R '@RG' reference_genome.fa fastq_r1 fastq_r2 > file.sam 
bwa_sam= bwa_storage_dir+r1_base_filename + ".sam"
print steps_to_process
if not os.path.isdir(bwa_storage_dir):
	subprocess.call(["mkdir", bwa_storage_dir])
if (not os.path.isfile(bwa_sam)) and("bwa-aligner" in steps_to_process):
	print "beginning alignment"
	print fastq_r1, fastq_r2
	cmd = BwaMem[:]
	cmd.append("-M")
	cmd.append("-t")
	cmd.append("6")
	cmd.append("-R")
	cmd.append(read_groups)
	#cmd.append(log)
	cmd.append(reference_genome)
	cmd.append(fastq_r1_location) 
	cmd.append(fastq_r2_location)
	print " ".join(cmd)
	f = open(bwa_sam, "w")
	errcode = subprocess.call(cmd, stdout=f)
	if(errcode == 0):
		print "BWA finished successfully"
	else:
		print "BWA failed !!!!"
		del steps_to_process[:]
	f.close()
	del cmd

# SAMFORMATCONVERTER
# =====================================================================
#SamFormatConverter I=bwa_output.sam O=converted_bwa_output.bamt
picard_bam = picard_storage_dir+r1_base_filename + ".bam"
if not os.path.isdir(picard_storage_dir):
	subprocess.call(["mkdir", picard_storage_dir])
if(not os.path.isfile(picard_bam)) and(("samformatconverter" in steps_to_process) or ("all" in steps_to_process)):
	cmd = SamFormatConverter[:]
	cmd.append("I=")
	cmd.append(bwa_sam)
	cmd.append("O=")
	cmd.append(picard_bam)
	cmd.append("MAX_RECORDS_IN_RAM=200000")
	cmd.append("TMP_DIR=")
	cmd.append(picard_tmp_dir)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "samformatconverter finished successfully"
	else:
		print "samformatconverter failed !!!!"
		del steps_to_process[:]
	os.remove(bwa_sam)
	del cmd


# SORTSAM
# =====================================================================
#SortSam I=input.bam O=sorted.bam SORT_ORDER=coordinate 
sortsam_bam = picard_storage_dir+r1_base_filename + "_sorted.bam"
if (not os.path.isfile(sortsam_bam)) and(("sortsam" in steps_to_process) or ("all" in steps_to_process)):
	print "sortsam_bam:",sortsam_bam
	cmd = SortSam[:]
	cmd.append("I=")
	cmd.append(picard_bam)
	cmd.append("O=")
	cmd.append(sortsam_bam)
	cmd.append("SORT_ORDER=coordinate")
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "sortsam finished successfully"
	else:
		print "sortsam failed !!!!"
		del steps_to_process[:]
	del cmd


# MARKDUPLICATES
# =====================================================================
#MarkDuplicates I=input.bam O=marked_duplicates.bam M=marked_dup_metrics.txt
marked_duplicates_bam= picard_storage_dir+"markduplicates/"+ r1_base_filename + "_marked_duplicates.bam"
marked_dup_metrics = picard_storage_dir+r1_base_filename+"_marked_dup_metrics.txt"
picard_tmp_dir = picard_storage_dir+"tmp/"
if ("markduplicates" in steps_to_process) or ("all" in steps_to_process):
	cmd = MarkDuplicates[:]
	cmd.append("I=")
	cmd.append(sortsam_bam)
	cmd.append("O=")
	cmd.append(marked_duplicates_bam)
	cmd.append("M=")
	cmd.append(marked_dup_metrics)
	cmd.append("TMP_DIR=")
	cmd.append(picard_tmp_dir)
        out_file = picard_storage_dir + r1_base_filename+"_markduplicates.log"
        err_file = picard_storage_dir + r1_base_filename+"_markduplicates.err"
        print " ".join(cmd)
        with open(err_file,"w")as outerr:
                with open(out_file,"w")as outf:
                        errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
                        if(errcode == 0):
                                print "markduplicates finished successfully"
                        else:
                                print "markduplicates failed !!!!"
                                del steps_to_process[:]
        del cmd

# MARKDUPLICATES (and remove)
# =====================================================================
#MarkDuplicates I=input.bam O=marked_duplicates.bam M=marked_dup_metrics.txt
remove_duplicates_bam= picard_storage_dir+"removeduplicates/"+ r1_base_filename + "_removed_duplicates.bam"
remove_dup_metrics = picard_storage_dir+"removeduplicates/"+r1_base_filename+"_removed_dup_metrics.txt"
picard_tmp_dir = picard_storage_dir+"tmp/"
if(not os.path.isfile(remove_duplicates_bam)) and (("removeduplicates" in steps_to_process) or ("all" in steps_to_process)):
	cmd = MarkDuplicates[:]
	cmd.append("I=")
	cmd.append(sortsam_bam)
	cmd.append("O=")
	cmd.append(remove_duplicates_bam)
	cmd.append("M=")
	cmd.append(remove_dup_metrics)
	cmd.append("TMP_DIR=")
	cmd.append(picard_tmp_dir)
	cmd.append("REMOVE_DUPLICATES=true")
        out_file = picard_storage_dir + r1_base_filename+"_markduplicates.log"
        err_file = picard_storage_dir + r1_base_filename+"_markduplicates.err"
        print " ".join(cmd)
        with open(err_file,"w")as outerr:
                with open(out_file,"w")as outf:
                        errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
                        if(errcode == 0):
                                print "markduplicates finished successfully"
                        else:
                                print "markduplicates failed !!!!"
                                del steps_to_process[:]
        del cmd

# COLLECTHSMETRICS
# =====================================================================
#CollectHsMetrics I=input.bam O=hsmetrics.txt R=reference_seq.fa BAIT_INTERVALS=bait_intervals.intervals TARGET_INTERVALS=target_intervals.intervals 
marked_duplicates_storage_bam= picard_storage_dir+"markduplicates/"+ r1_base_filename + "_marked_duplicates.bam"
hs_metrics= alignment_qc_storage_dir+"collecthsmetrics/"+r1_base_filename + "_hsmetrics.txt"
bait_intervals = "/dataVolume/storage/rb_pipeline/agilent_baits_hg19.interval_list"
target_intervals = "/dataVolume/storage/rb_pipeline/agilent_targets_hg19.interval_list"
if(not os.path.isfile(hs_metrics)) and (("collecthsmetrics" in steps_to_process) or ("all" in steps_to_process)):
	cmd = CollectHsMetrics[:]
	cmd.append("I=")
	cmd.append(marked_duplicates_storage_bam)
	cmd.append("O=") 
	cmd.append(hs_metrics)
	cmd.append("R=")
	cmd.append(reference_genome)
	cmd.append("BAIT_INTERVALS=")
	cmd.append(bait_intervals)
	cmd.append("TARGET_INTERVALS=")
	cmd.append(target_intervals)
	out_file = alignment_qc_storage_dir + r1_base_filename+"_hsmetrics.log"
	err_file = alignment_qc_storage_dir + r1_base_filename+"_hsmetrics.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
		 	errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
			if(errcode == 0):
				print "collecthsmetrics finished successfully"
			else:
				print "collecthsmetrics failed !!!!"
				del steps_to_process[:]
	del cmd
	
# COLLECTINSERTSIZEMETRICS
# =====================================================================
#CollectInsertSizeMetrics I=input.bam O=insertsizemetrics.txt R=reference_seq.fa BAIT_INTERVALS=bait_intervals.intervals TARGET_INTERVALS=target_intervals.intervals 
is_metrics= alignment_qc_storage_dir+"collectinsertsizemetrics/"+r1_base_filename + "_ismetrics.txt"
is_metrics_histogram= alignment_qc_storage_dir+"collectinsertsizemetrics/"+r1_base_filename + "_ismetrics_histogram.pdf"

if(not os.path.isfile(is_metrics)) and (("collectismetrics" in steps_to_process) or ("alignmentqc" in steps_to_process)):
	cmd = CollectIsMetrics[:]
	cmd.append("I=")
	cmd.append(marked_duplicates_storage_bam)
	cmd.append("O=") 
	cmd.append(is_metrics)
	cmd.append("H=")
	cmd.append(is_metrics_histogram)
	out_file = alignment_qc_storage_dir + r1_base_filename+"_ismetrics.log"
	err_file = alignment_qc_storage_dir + r1_base_filename+"_ismetrics.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
		 	errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
			if(errcode == 0):
				print "collectismetrics finished successfully"
			else:
				print "collectismetrics failed !!!!"
				del steps_to_process[:]
	del cmd
	
# COLLECTALIGNMENTSUMMARYMETRICS
# =====================================================================
#CollectAlignmentSummaryMetrics I=input.bam O=alignsummetrics.txt R=reference_seq.fa 
alignsum_metrics= alignment_qc_storage_dir+"collectalignmentsummarymetrics/"+r1_base_filename + "_alignsummetrics.txt"
if(not os.path.isfile(alignsum_metrics)) and (("collectalignmentsummarymetrics" in steps_to_process) or ("alignmentqc" in steps_to_process)):
	cmd = CollectAlignmentSummaryMetrics[:]
	cmd.append("R=")
	cmd.append(reference_genome)
	cmd.append("I=")
	cmd.append(marked_duplicates_storage_bam)
	cmd.append("O=") 
	cmd.append(alignsum_metrics)
	out_file = alignment_qc_storage_dir + r1_base_filename+"_alignment_summary_metrics.log"
	err_file = alignment_qc_storage_dir + r1_base_filename+"_alignment_summary_metrics.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
		 	errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
			if(errcode == 0):
				print "collectalignmentsummarymetrics finished successfully"
			else:
				print "collectalignmentsummarymetrics failed !!!!"
				del steps_to_process[:]
	del cmd



# BUILDBAMINDEX
# =====================================================================
#BuildBamIndex I=input.bam 
bam_index= picard_storage_dir+"markduplicates/"+r1_base_filename+ "_marked_duplicates.bai"
if ("buildbamindex" in steps_to_process):
	cmd = BuildBamIndex[:]
	cmd.append("I=")
	cmd.append(marked_duplicates_bam)
	cmd.append("O=")
	cmd.append(bam_index)
	print " ".join(cmd)
	out_file = picard_storage_dir+"markduplicates/" + r1_base_filename+"_buildbamindex.log"
	err_file = picard_storage_dir+"markduplicates/" + r1_base_filename+"_buildbamindex.err"
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
                	errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
			if(errcode == 0):
                        	print "buildbamindex finished successfully"
                        else:
                                print "buildbamindex failed !!!!"
                                del steps_to_process[:]
        del cmd

"""
#RealignerTargetCreator and IndelRealigner are not necessary if running Mutect2 or haplotypecaller!

# REALIGNERTARGETCREATOR
# =====================================================================
#RealignerTargetCreator -T=RealignerTargetCreator -R=reference.fasta -I=input.bam --known=indels.vcf -o=forIndelRealigner.intervals
intervals_for_IR = gatk_storage_dir+r1_base_filename+"_forIndelRealigner.intervals"
indels_vcf_Mills = "/dataVolume/storage/rb_pipeline/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
indels_vcf_1000g = "/dataVolume/storage/rb_pipeline/1000G_phase1.indels.hg19.sites.vcf"
if not os.path.isdir(gatk_storage_dir):
	subprocess.call(["mkdir", gatk_storage_dir])
if(not os.path.isfile(intervals_for_IR)) and (("realignertargetcreator" in steps_to_process) or ("all" in steps_to_process)):
	cmd = RealignerTargetCreator[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I")
	cmd.append(marked_duplicates_bam)
	cmd.append("-known")
	cmd.append(indels_vcf_Mills)
	cmd.append("-known")
	cmd.append(indels_vcf_1000g)
	cmd.append("-o")
	cmd.append(intervals_for_IR)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "realignertargetcreator finished successfully"
	else:
		print "realignertargetcreator failed !!!!"
		del steps_to_process[:]
	del cmd
	
# INDELREALIGNER
# =====================================================================
#IndelRealigner -T=IndelRealigner -R=reference.fasta -I=input.bam -known=indels.vcf -targetIntervals=intervalListFromRTC.intervals -o=realignedBam.bam
realigned_bam = gatk_storage_dir+r1_base_filename+"_realigned.bam"
if(not os.path.isfile(realigned_bam)) and (("indelrealigner" in steps_to_process) or ("all" in steps_to_process)):
	cmd = IndelRealigner[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I")
	cmd.append(marked_duplicates_bam)
	cmd.append("-known")
	cmd.append(indels_vcf_Mills)
	cmd.append("-known")
	cmd.append(indels"_vcf_1000g)
	cmd.append("-targetIntervals")
	cmd.append(intervals_for_IR)
	cmd.append("-o")
	cmd.append(realigned_bam)
	print " ".join(cmd)
	errcode = subprmocess.call(cmd)
	if(errcode == 0):
		print "indelrealigner finished successfully"
	else:
		print "indelrealigner failed !!!!"
		del steps_to_process[:] 
	del cmd
"""

# BASERECALIBRATOR
# =====================================================================
#BaseRecalibrator -T=IndelRealigner -R=reference.fasta -I=input.bam -known=indels.vcf -targetIntervals=intervalListFromRTC.intervals -o=realignedBam.bam 
realigned_bam = gatk_storage_dir+r1_base_filename+"_realigned.bam"
recal_report = gatk_storage_dir+r1_base_filename+"_recal_report.table"
indels_vcf_dbsnp = "/dataVolume/storage/rb_pipeline/dbsnp_138.hg19.excluding_sites_after_129.vcf"
indels_vcf_Mills = "/dataVolume/storage/rb_pipeline/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
indels_vcf_1000g = "/dataVolume/storage/rb_pipeline/1000G_phase1.indels.hg19.sites.vcf"
if not (os.path.isfile(recal_report)) and (("baserecalibrator" in steps_to_process) or ("all" in steps_to_process)):
	cmd = BaseRecalibrator[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I")
	cmd.append(marked_duplicates_bam)
	cmd.append("-knownSites")
	cmd.append(indels_vcf_Mills)
	cmd.append("-knownSites")
	cmd.append(indels_vcf_1000g)
	cmd.append("-knownSites")
	cmd.append(indels_vcf_dbsnp)
	cmd.append("-o")
	cmd.append(recal_report)
	cmd.append("-rf")
	cmd.append("BadCigar")
	out_file = gatk_storage_dir + r1_base_filename+"_baserecalibrator.log"
        err_file = gatk_storage_dir + r1_base_filename+"_baserecalibrator.err"
        print " ".join(cmd)
        with open(err_file,"w")as outerr:
                with open(out_file,"w")as outf:
                        errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
                        if(errcode == 0):
                                print "baserecalibrator finished successfully"
                        else:
                                print "baserecalibrator failed !!!!"
                                del steps_to_process[:]
        del cmd

"""
# ANALYZECOVARIATES
# =====================================================================
#AnalyzeCovariates -T=AnalyzeCovariates -R=reference.fasta -I=input.bam -BQSR= recalibration_report.grp -o=output.bam
recal_plot = gatk_storage_dir+r1_base_filename+"_recalQC.pdf"
if(not os.path.isfile(recal_plot)) and (("analyzecovariates" in steps_to_process) or ("all" in steps_to_process)):
	cmd = AnalyzeCovariates[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-BQSR")
	cmd.append(recal_report)
	cmd.append("-plots")
	cmd.append(recal_plot)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "analyzecovariates finished successfully"
	else:
		print "analyzecovariates failed !!!!"
		del steps_to_process[:]
	del cmd
"""
	
# PRINTREADS
# =====================================================================
#PrintReads -T=PrintReads -R=reference.fasta -I=input.bam -BQSR= recalibration_report.grp -o=output.bam
recal_bam = gatk_storage_dir+r1_base_filename+"_recalibrated.bam"
if (not os.path.isfile(recal_bam)) and (("printreads" in steps_to_process) or ("all" in steps_to_process)):
	cmd = PrintReads[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I")
	cmd.append(marked_duplicates_bam)
	cmd.append("-BQSR")
	cmd.append(recal_report)
	cmd.append("-o")
	cmd.append(recal_bam)
	out_file = gatk_storage_dir + r1_base_filename+"_printreads.log"
	err_file = gatk_storage_dir + r1_base_filename+"_printreads.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
		 	errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
			if(errcode == 0):
				print "printreads finished successfully"
			else:
				print "printreads failed !!!!"
				del steps_to_process[:]
	del cmd
	
# recalibrated BUILDBAMINDEX
# =====================================================================
#BuildBamIndex I=input.bam 
recal_bam_index= gatk_storage_dir+r1_base_filename+ "_recalibrated.bai"
if ("recalibratedbamindex" in steps_to_process):
	cmd = BuildBamIndex[:]
	cmd.append("I=")
	cmd.append(recal_bam)
	cmd.append("O=")
	cmd.append(recal_bam_index)
	out_file = gatk_storage_dir + r1_base_filename+"_recalibratedindex.log"
        err_file = gatk_storage_dir + r1_base_filename+"_recalibratedindex.err"
        print " ".join(cmd)
        with open(err_file,"w")as outerr:
                with open(out_file,"w")as outf:
                        errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
                        if(errcode == 0):
                                print "recalibratedbamindex finished successfully"
                        else:
                                print "recalibratedbamindex failed !!!!"
                                del steps_to_process[:]
        del cmd

# MUTECT2
# =====================================================================
mutect2_vcf=mutect2_storage_dir+r1_base_filename+"_mutect2.vcf"
if not os.path.isdir(mutect2_storage_dir):
	subprocess.call(["mkdir", mutect2_storage_dir])
if "-T" in recal_bam:
	primary_recal_bam = recal_bam
	match_normal_recal_bam = re.sub('-T_1', '-N_1', primary_recal_bam)
elif "-CL" in recal_bam:
	primary_recal_bam = recal_bam
	match_normal_recal_bam = re.sub('-CL_1', '-N_1', primary_recal_bam)
else:
	steps_to_process += "skip_mutect2"
if(not os.path.isfile(mutect2_vcf)) and ("skip_mutect2" not in steps_to_process) and(("mutect2" in steps_to_process) or ("all" in steps_to_process)):
	cmd  = Mutect2[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I:tumor") # tumor or cell line depending on input sample!
	cmd.append(primary_recal_bam)
	cmd.append("-I:normal")
	cmd.append(match_normal_recal_bam)
	cmd.append("--dbsnp")
	cmd.append(indels_vcf_dbsnp)
	cmd.append("-L")
	cmd.append(target_intervals) #need to evaluate 
	cmd.append("-o")
	cmd.append(mutect2_vcf)
	out_file = gatk_storage_dir + r1_base_filename+"_mutect2.log"
	err_file = gatk_storage_dir + r1_base_filename+"_mutect2.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
		 	errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
			if(errcode == 0):
				print "mutect2 finished successfully"
			else:
				print "mutect2 failed !!!!"
				del steps_to_process[:]
	del cmd

# MUTECT2-PANEL OF NORMALS
# =====================================================================
pon_vcf=mutect2_storage_dir+r1_base_filename+"_pon_mutect2.vcf"
if not os.path.isdir(mutect2_storage_dir):
        subprocess.call(["mkdir", mutect2_storage_dir])
if(not os.path.isfile(mutect2_vcf)) and(("panelofnormals" in steps_to_process) or ("all" in steps_to_process)):
        cmd  = Mutect2[:]
        cmd.append("-R")
        cmd.append(reference_genome)
        cmd.append("-I:tumor") # tumor or cell line depending on input sample!
        cmd.append(recal_bam)
        cmd.append("--artifact_detection_mode")
        cmd.append("-L")
        cmd.append(target_intervals) #need to evaluate 
        cmd.append("-o")
        cmd.append(pon_vcf)
        out_file = storage_directory + r1_base_filename+"_pon.log"
        err_file = storage_directory + r1_base_filename+"_pon.err"
        print " ".join(cmd)
        with open(err_file,"w")as outerr:
                with open(out_file,"w")as outf:
                        errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
                        if(errcode == 0):
                                print "panel of normals finished successfully"
                        else:
                                print "panel of normals failed !!!!"
                                del steps_to_process[:]
        del cmd

# MUTECT2SWAP	
# =====================================================================
mutect2_swap_vcf=mutect2_storage_dir+r1_base_filename+"_mutect2_swap.vcf"
if not os.path.isdir(mutect2_storage_dir):
	subprocess.call(["mkdir", mutect2_storage_dir])
if "-T" in recal_bam:
	match_normal_recal_bam = recal_bam
	primary_recal_bam = re.sub('-T_1', '-N_1', match_normal_recal_bam)
elif "-CL" in recal_bam:
	match_normal_recal_bam = recal_bam
	primary_recal_bam = re.sub('-CL_1', '-N_1', match_normal_recal_bam)
else:
	steps_to_process += "skip_mutect2"
if(not os.path.isfile(mutect2_swap_vcf)) and ("skip_mutect2" not in steps_to_process) and(("mutect2swap" in steps_to_process) or ("all" in steps_to_process)):
	cmd  = Mutect2[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I:tumor") # swapped tumor\normal in this command!
	cmd.append(primary_recal_bam)
	cmd.append("-I:normal")
	cmd.append(match_normal_recal_bam)
	cmd.append("--dbsnp")
	cmd.append(indels_vcf_dbsnp)
	cmd.append("-L")
	cmd.append(target_intervals) #need to evaluate 
	cmd.append("-o")
	cmd.append(mutect2_swap_vcf)
	out_file = mutect2_storage_dir + r1_base_filename+"_mutect2_swap.log"
	err_file = mutect2_storage_dir + r1_base_filename+"_mutect2_swap.err"
	print " ".join(cmd)
	with open(err_file,"w")as outerr:
		with open(out_file,"w")as outf:
		 	errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
			if(errcode == 0):
				print "mutect2swap finished successfully"
			else:
				print "mutect2swap failed !!!!"
				del steps_to_process[:]
	del cmd

# MUTECT (for troubleshooting low SNP count in cobrinik data
# =====================================================================
mutect_vcf=mutect_storage_dir+r1_base_filename+"_mutect.vcf"
mutect_coverage=mutect_storage_dir+r1_base_filename+"_mutect_coverage.wig.txt"
regions_bed="/dataVolume/storage/rb_pipeline/Agilent_Sureselect_V5_Regions.bed"
if not os.path.isdir(mutect_storage_dir):
	subprocess.call(["mkdir", mutect_storage_dir])
if "-T" in recal_bam:
	primary_recal_bam = recal_bam
	match_normal_recal_bam = re.sub('-T_1', '-N_1', primary_recal_bam)
elif "-CL" in recal_bam:
	primary_recal_bam = recal_bam
	match_normal_recal_bam = re.sub('-CL_1', '-N_1', primary_recal_bam)
else:
	steps_to_process += "skip_mutect"
if(not os.path.isfile(mutect_vcf)) and ("skip_mutect" not in steps_to_process) and("mutect" in steps_to_process):
	cmd  = Mutect[:]
	cmd.append("--analysis_type")
	cmd.append("MuTect")
	cmd.append("--reference_sequence")
	cmd.append(reference_genome)
	cmd.append("--dbsnp")
	cmd.append(indels_vcf_dbsnp)
	cmd.append("--intervals")
	cmd.append(target_intervals)
	cmd.append("--input_file:normal")
	cmd.append(match_normal_recal_bam)
	cmd.append("--input_file:tumor")
	cmd.append(primary_recal_bam)
	cmd.append("--out")
	cmd.append(mutect_vcf)
	cmd.append("--coverage_file")
	cmd.append(mutect_coverage)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "mutect finished successfully"
	else:
		print "mutect failed !!!!"
		del steps_to_process[:]
	del cmd

# HAPLOTYPECALLER
# =====================================================================
if not os.path.isdir(haplocaller_storage_dir):
	subprocess.call(["mkdir", haplocaller_storage_dir])
haplocaller_vcf = haplocaller_storage_dir+r1_base_filename+"_haplocaller.vcf"
if (not os.path.isfile(haplocaller_vcf)) and ("skip haplotypecaller" not in steps_to_process) and ("haplocaller" in steps_to_process):
	cmd  = HaplotypeCaller[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I")
	cmd.append(primary_recal_bam)
	cmd.append("--dbsnp")
	cmd.append(indels_vcf_dbsnp)
	cmd.append("-stand_call_conf")
	cmd.append("30")
	cmd.append("-L")
	cmd.append(target_intervals)
	cmd.append("-o")
	cmd.append(haplocaller_vcf)
        out_file = gatk_storage_dir + r1_base_filename+"_haplocaller.log"
        err_file = gatk_storage_dir + r1_base_filename+"_haplocaller.err"
        print " ".join(cmd)
        with open(err_file,"w")as outerr:
                with open(out_file,"w")as outf:
                        errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
                        if(errcode == 0):
                                print "haplocaller finished successfully"
                        else:
                                print "haplocaller failed !!!!"
                                del steps_to_process[:]
        del cmd


# VARIANTFILTRATION
# =====================================================================
#VariantFiltration -T=VariantFiltration -R=reference.fasta -o=output.vcf --variant= input.vcf --filterExpression= "AB < 0.2 || MQ0 > 50" --filterName= "SomeFilterName"
RB_vcf = mutect2_storage_dir+r1_base_filename+"_filtered.vcf"
if(not os.path.isfile(RB_vcf)) and (("variantfiltration" in steps_to_process) or ("all" in steps_to_process)):
	cmd = VariantFiltration[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-o")
	cmd.append(RB_vcf)
	cmd.append("--variant")
	cmd.append(mutect2_vcf)
	cmd.append("--clusterWindowSize") 
	cmd.append("10")
	cmd.append("--filterExpression")
	cmd.append("MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)")
	cmd.append("--filterName") 
	cmd.append("HARD_TO_VALIDATE")
	cmd.append("--filterExpression") 
	cmd.append("DP < 5 ")
	cmd.append("--filterName") 
	cmd.append("LowCoverage")
	cmd.append("--filterExpression") 
	cmd.append("QUAL < 30.0 ")
	cmd.append("--filterName")
	cmd.append("VeryLowQual")
	cmd.append("--filterExpression") 
	cmd.append("QUAL >= 30.0 && QUAL < 50.0 ")
	cmd.append("--filterName") 
	cmd.append("LowQual")
	cmd.append("--filterExpression")
	cmd.append("QD < 1.5 ")
	cmd.append("--filterName") 
	cmd.append("LowQD")
	cmd.append("--filterExpression") 
	cmd.append("FS > 60.0 ")
	cmd.append("--filterName") 
	cmd.append("FisherStrandBiasSNP")
	cmd.append("--filterExpression") 
	cmd.append("FS > 200.0 ")
	cmd.append("--filterName") 
	cmd.append("FisherStrandBiasINDEL")
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "variantfiltration finished successfully"
	else:
		print "variantfiltration failed !!!!"
		del steps_to_process[:]
	del cmd
	
	# SELECTVARIANTS
# =====================================================================
#SelectVariants -T=VariantFiltration -R=reference.fasta -o=output.vcf --variant= input.vcf --filterExpression= "AB < 0.2 || MQ0 > 50" --filterName= "SomeFilterName"
selected_vcf = mutect2_storage_dir+r1_base_filename+"_selected.vcf"
if(not os.path.isfile(selected_vcf)) and (("selectvariants" in steps_to_process) or ("all" in steps_to_process)):
	cmd = SelectVariants[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("--variant")
	cmd.append(RB_vcf)
	cmd.append("-o")
	cmd.append(selected_vcf)
	cmd.append("-select")
	cmd.append("vc.isNotFiltered()")
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "selectvariants finished successfully"
	else:
		print "selectvariants failed !!!!"
		del steps_to_process[:]
	del cmd
	
	
# ANNOVAR
# =====================================================================
RB_avinput = annovar_storage_dir+r1_base_filename+".avinput"
RB_variants = annovar_storage_dir
#convert vcf file to annovar input file format
if not os.path.isdir(annovar_storage_dir):
	subprocess.call(["mkdir", annovar_storage_dir])
if(not os.path.isfile(RB_avinput)) and (("annovar" in steps_to_process) or ("all" in steps_to_process)):
	cmd1 = ConvertToAnnovar[:]
	cmd1.append("-format")
	cmd1.append("vcf4")
	cmd1.append(RB_vcf)
	cmd1.append("-outfile")
	cmd1.append(RB_avinput)
	cmd1.append("-includeinfo")
#	cmd1.append("-withzyg")
	cmd1.append("-withfreq")
	cmd1.append("-include")
	cmd1.append("-comment")
	print " ".join(cmd1)
	errcode1 = subprocess.call(cmd1)
	if(errcode1 == 0):
		print "annovar finished successfully"
	else:
		print "annovar failed !!!!"
		del steps_to_process[:]
	del cmd1

	cmd2 = GeneAnnotation[:]
	cmd2.append("-out")
	cmd2.append(r1_base_filename)
	cmd2.append("-build")
	cmd2.append("hg19")
	cmd2.append(RB_avinput)
	cmd2.append("/dataVolume/storage/rb_pipeline/humandb/")
	print " ".join(cmd2)
	errcode2 = subprocess.call(cmd2)
	if(errcode2 == 0):
		print "annovar finished successfully"
	else:
		print "annovar failed !!!!"
		del steps_to_process[:]
	del cmd2
	# need also region and filter-based annotations?

	cmd3 = FilterAnnotation[:]
	cmd3.append(RB_avinput)
	cmd3.append("/dataVolume/storage/rb_pipeline/humandb/")
	cmd3.append("-protocol")
	cmd3.append("dbnsfp30a")
	cmd3.append("-operation")
	cmd3.append("f")
	cmd3.append("-build")
	cmd3.append("hg19")
	cmd3.append("-nastring")
	print " ".join(cmd3)
	errcode3 = subprocess.call(cmd3)
	if(errcode3 == 0):
			print "annovar filterannotation finished successfully"
	else:
			print "annovar filterannotation failed !!!!"
			del steps_to_process[:]
	del cmd3


	#how to include information from pathogen predicting programs?

# COPYWRITER
# =====================================================================
copywriter_input_dir = "/dataVolume/storage/rb_pipeline/picard/removeduplicates/"+r1_base_filename
trio_name = r1_base_filename.rsplit('-', 1)[0]
copywriter_output_dir = output_directory+"copywriter/"+r1_base_filename
if not os.path.isdir(copywriter_output_dir):
          subprocess.call(["mkdir", copywriter_output_dir])     
if ("copywriter" in steps_to_process) or ("all" in steps_to_process):
        cmd = CopyWriter[:]
        cmd.append(copywriter_input_dir)
        cmd.append(trio_name)
        cmd.append(r1_base_filename)
        cmd.append(copywriter_output_dir)
	cmd.append("50")
#        out_file = output_directory+r1_base_filename+ "_copywriter.log"
#        err_file = output_directory+r1_base_filename+"_copywriter.err"
        print " ".join(cmd)
#        with open(err_file,"w")as outerr:
#                with open(out_file,"w")as outf:
#                        errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
	errcode = subprocess.call(cmd)
        if(errcode == 0):
                print "copywriter finished successfully"
        else:
                print "copywriter failed !!!!"
                del steps_to_process[:]
        del cmd

# COPYWRITER (HI RES)
# =====================================================================
copywriter_input_dir = "/dataVolume/storage/rb_pipeline/picard/removeduplicates/"+r1_base_filename
trio_name = r1_base_filename.rsplit('-', 1)[0]
hi_res_copywriter_output_dir = output_directory+"hi_res_copywriter/"+r1_base_filename
if not os.path.isdir(hi_res_copywriter_output_dir):
          subprocess.call(["mkdir", hi_res_copywriter_output_dir])
print hi_res_copywriter_output_dir     
if ("hirescopywriter" in steps_to_process) or ("all" in steps_to_process):
        cmd = HiResCopyWriter[:]
        cmd.append(copywriter_input_dir)
        cmd.append(trio_name)
        cmd.append(r1_base_filename)
        cmd.append(hi_res_copywriter_output_dir)
#        out_file = output_directory+r1_base_filename+ "_hi_res_copywriter.log"
#        err_file = output_directory+r1_base_filename+"_hi_res_copywriter.err"
        print " ".join(cmd)
#        with open(err_file,"w")as outerr:
#                with open(out_file,"w")as outf:
#                        errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
	errcode = subprocess.call(cmd)
        if(errcode == 0):
                print "hi-res-copywriter finished successfully"
        else:
                print "hi-res-copywriter failed !!!!"
                del steps_to_process[:]
        del cmd	

#VCFTOOLS
# =====================================================================
#if(not os.path.isfile(RB_vcf)) and (("vcftools" in steps_to_process) or ("all" in steps_to_process)):
#cmd = ["vcftools"]
#cmd.append("--vcf")
#cmd.append(RB_vcf)
#cmd.append("--out")
#cmd.append
#	del steps_to_process[:]
#	del cmd

#SAMTOOLS MPILEUP
# =====================================================================
#if(not os.path.isfile(RB_vcf)) and (("vcftools" in steps_to_process) or ("all" in steps_to_process)):

#BAFSEGEMENTATION
# =====================================================================
#if(not os.path.isfile(RB_vcf)) and (("vcftools" in steps_to_process) or ("all" in steps_to_process)):







