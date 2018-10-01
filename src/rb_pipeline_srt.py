#!/usr/bin/python

# argv[1,2] = R1/R2.fastq files
# argv[3] = reference genome (*.fasta)
# argv[4] = step to process

import subprocess
import sys
import re
import os
import gzip 


num_threads = "6"

# *********************************************************************
# DEFINITION OF PATHS

fastq_r1_location = sys.argv[1]
fastq_r2_location = sys.argv[2]
reference_genome = sys.argv[3]
steps_to_process = sys.argv[4].split(",")

Fastqc = ["fastqc"]
BwaMem = ["bwa", "mem"]
SamFormatConverter = ["java",  "-jar", "/usr/share/java/picard.jar", "SamFormatConverter"]
SortSam  = ["java", "-jar", "/usr/share/java/picard.jar", "SortSam"]
MarkDuplicates = ["java", "-jar", "/usr/share/java/picard.jar", "MarkDuplicates"]
CollectHsMetrics = ["java", "-jar", "/usr/share/java/picard.jar", "CollectHsMetrics"]
BuildBamIndex = ["java", "-jar", "/usr/share/java/picard.jar", "BuildBamIndex"]
RealignerTargetCreator = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "RealignerTargetCreator"]
IndelRealigner = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "IndelRealigner"]
BaseRecalibrator = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "BaseRecalibrator"]
PrintReads = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "PrintReads"]
AnalyzeCovariates = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "AnalyzeCovariates"]
HaplotypeCaller = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "HaplotypeCaller"]
GenotypeGVCFs = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "GenotypeGVCFs"]
VariantFiltration = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK.jar", "-T", "VariantFiltration"]
Trimmomatic = ["java", "-jar", "/usr/share/java/trimmomatic.jar"]
ConvertToAnnovar = ["/usr/share/perl/5.22.1/annovar/convert2annovar.pl"]
GeneAnnotation = ["/usr/share/perl/5.22.1/annovar/annotate_variation.pl"]
CopyWriter = ["Rscript","/media/sf_RB_exome_project/CopywriteR_exomes/CopywriteR_test_script.R"]
illuminaclip = "TruSeq_DNA_v1_v1_LT_multiplex_allInOneIndex1_22.fa:2:30:10"
min_base_quality = "3"
sliding_window = "4:15"
min_read_length  = "36"

output_directory = "/dataVolume/storage/rb_pipeline/"
picard_output_dir = output_directory+"picard/"
bwa_output_dir = output_directory+"bwa/"
picard_tmp_dir = output_directory+"tmp/"
gatk_output_dir = output_directory + "gatk/"
haplocaller_output_dir = output_directory+ "gatk/haplocaller/"
annovar_output_dir = output_directory+"annovar/"
bafsegment_output_dir = output_directory+"bafseqgment/"
r1_base_filename = fastq_r1_location.replace(".fastq.gz","").split("/")[-1]
r2_base_filename = fastq_r2_location.replace(".fastq.gz","").split("/")[-1]
fastq_r1 = fastq_r1_location.split("/")[-1]
fastq_r2 = fastq_r2_location.split("/")[-1]
RBID = r1_base_filename.split("_")[0]

# /DEFINITION OF PATHS
# *********************************************************************

# READ GROUP PARSING
# =====================================================================

with gzip.open(fastq_r1_location) as fasta:
	s = fasta.readline().split(":")
	RGID =  "_".join(s[2:5])
	RGSM = r1_base_filename
	RGPL = "Illumina"
	read_groups = "@RG\\tID:"+RGID+"\\tSM:"+RGSM+"\\tPL:"+RGPL
	print read_groups
	
# TRIMMOMATIC
# =====================================================================
trimmomatic_output_dir = output_directory + "trimmomatic/"

fastq_r1_trimmed = trimmomatic_output_dir+ r1_base_filename + ".trimmed.fastq"
fastq_r2_trimmed = trimmomatic_output_dir+ r2_base_filename + ".trimmed.fastq"
fastq_r1_trimmed_unpaired = trimmomatic_output_dir+ r1_base_filename + ".trimmed.unpaired.fastq"
fastq_r2_trimmed_unpaired = trimmomatic_output_dir+ r2_base_filename + ".trimmed.unpaired.fastq"
if(not os.path.isfile(fastq_r1_trimmed)) and(("trimmomatic" in steps_to_process) or ("all" in steps_to_process)):
	print sys.argv[4]
	subprocess.call(["mkdir",trimmomatic_output_dir])
	print "running trimmomatic"
	cmd = Trimmomatic[:]
	#print cmd
	cmd.append("PE")
	cmd.append("-threads")
	cmd.append("4")
	#cmd.append("-trimlog")
	#cmd.append(log)
	cmd.append(fastq_r1_location)
	cmd.append(fastq_r2_location)
	cmd.append(fastq_r1_trimmed)
	cmd.append(fastq_r1_trimmed_unpaired)
	cmd.append(fastq_r2_trimmed)
	cmd.append(fastq_r2_trimmed_unpaired)
	cmd.append("ILLUMINACLIP:"+illuminaclip)
	cmd.append("LEADING:"+min_base_quality)
	cmd.append("TRAILING:"+min_base_quality)
	cmd.append("SLIDINGWINDOW:"+sliding_window)
	cmd.append("MINLEN:"+min_read_length)
	print " ".join(cmd)
	#trimmomatic = subprocess.Popen(cmd)
	#trimmomatic.wait()
	#print trimmomatic.returncode
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "trimmomatic for sample "+r1_base_filename+" finished successfully"
	else:
		print "fastqc for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd

# FASTQC
# =====================================================================
fastqc_output_dir = output_directory + "fastqc/"
fastqc_output = fastqc_output_dir+r1_base_filename+"_fasqc.zip"
if(not os.path.isfile(fastqc_output)) and(("fastqc" in steps_to_process) or ("all" in steps_to_process)):
	subprocess.call(["mkdir", fastqc_output_dir])
	print "running fastqc on the trimmed file"
	cmd = Fastqc[:]
	cmd.append(fastq_r1_trimmed)
	cmd.append(fastq_r2_trimmed)
	cmd.append(fastq_r1_location)
	cmd.append(fastq_r2_location)
	cmd.append("--outdir="+fastqc_output_dir)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "fastqc for sample "+r1_base_filename+" finished successfully"
	else:
		print "fastqc for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd
	subprocess.call(["mkdir", bwa_output_dir])

# BWA ALIGNER
# =====================================================================
#BWA MEM -aM -t 6 -R '@RG' reference_genome.fa fastq_r1 fastq_r2 > file.sam 
bwa_sam= bwa_output_dir+r1_base_filename+ ".trimmed" + ".sam"
if (not os.path.isfile(bwa_sam)) and ("bwa-aligner" in steps_to_process):
	print "beginning alignment"
	print fastq_r1, fastq_r2
	cmd = BwaMem[:]
	cmd.append("-M")
	cmd.append("-t")
	cmd.append("4")
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
		print "BWA for sample "+r1_base_filename+" finished successfully"
	else:
		print "BWA for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	f.close()
	del cmd
	subprocess.call(["mkdir", picard_output_dir])

# SAMFORMATCONVERTER
# =====================================================================
#SamFormatConverter I=bwa_output.sam O=converted_bwa_output.bamt
picard_bam = picard_output_dir+r1_base_filename + ".bam"
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
		print "samformatconverter for sample "+r1_base_filename+" finished successfully"
	else:
		print "samformatconverter for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd
	os.remove(bwa_sam)


# SORTSAM
# =====================================================================
#SortSam I=input.bam O=sorted.bam SORT_ORDER=coordinate 
sortsam_bam = picard_output_dir+ r1_base_filename + "_sorted.bam"
if (not os.path.isfile(sortsam_bam)) and(("sortsam" in steps_to_process) or ("all" in steps_to_process)):
	subprocess.call(["mkdir",picard_output_dir])
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
		print "sortsam for sample "+r1_base_filename+" finished successfully"
	else:
		print "sortsam for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd

# MARKDUPLICATES
# =====================================================================
#MarkDuplicates I=input.bam O=marked_duplicates.bam M=marked_dup_metrics.txt
marked_duplicates_bam= picard_output_dir+"markduplicates/"+ r1_base_filename + "_marked_duplicates.bam"
marked_dup_metrics = picard_output_dir+"markduplicates/"+r1_base_filename+"_marked_dup_metrics.txt"
picard_tmp_dir = picard_output_dir+"tmp/"
if(not os.path.isfile(marked_duplicates_bam)) and (("markduplicates" in steps_to_process) or ("all" in steps_to_process)):
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
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "markduplicates for sample "+r1_base_filename+" finished successfully"
	else:
		print "markduplicates for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd


# COLLECTHSMETRICS
# =====================================================================
#CollectHsMetrics I=input.bam O=hsmetrics.txt R=reference_seq.fa BAIT_INTERVALS=bait_intervals.intervals TARGET_INTERVALS=target_intervals.intervals 
hs_metrics= picard_output_dir+"collecthsmetrics/"+r1_base_filename + "_hsmetrics.txt"
bait_intervals = "SeqCapEZ_Exomev3-nonOverlapping-PT.intervals"
target_intervals = "refGene-21Nov2014-CDSnonoverlapping-PT.intervals"
if(not os.path.isfile(hs_metrics)) and (("collecthsmetrics" in steps_to_process) or ("all" in steps_to_process)):
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
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "collecthsmetrics for sample "+r1_base_filename+" finished successfully"
	else:
		print "collecthsmetrics for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd
	
# COLLECTINSERTSIZEMETRICS
# =====================================================================
#CollectInsertSizeMetrics I=input.bam O=insertsizemetrics.txt R=reference_seq.fa BAIT_INTERVALS=bait_intervals.intervals TARGET_INTERVALS=target_intervals.intervals 
is_metrics= picard_output_dir+"collectinsertsizemetrics/"+r1_base_filename + "_ismetrics.txt"
if(not os.path.isfile(is_metrics)) and (("collectismetrics" in steps_to_process) or ("all" in steps_to_process)):
	cmd = CollectIsMetrics[:]
	cmd.append("I=")
	cmd.append(marked_duplicates_bam)
	cmd.append("O=") 
	cmd.append(is_metrics)
	cmd.append("R=")
	cmd.append(reference_genome)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "collectismetrics for sample "+r1_base_filename+" finished successfully"
	else:
		print "collectismetrics for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd
	
# COLLECTALIGNMENTSUMMARYMETRICS
# =====================================================================
#CollectAlignmentSummaryMetrics I=input.bam O=alignsummetrics.txt R=reference_seq.fa 
alignsum_metrics= picard_output_dir+"alignment_qc/collectalignmentsummarymetrics/"+r1_base_filename + "_alignsummetrics.txt"
if(not os.path.isfile(alignsum_metrics)) and (("collectalignmentsummarymetrics" in steps_to_process) or ("all" in steps_to_process)):
	cmd = CollectAlignmentSummaryMetrics[:]
	cmd.append("I=")
	cmd.append(marked_duplicates_bam)
	cmd.append("O=") 
	cmd.append(alignsum_metrics)
	cmd.append("R=")
	cmd.append(reference_genome)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "collectalignmentsummarymetrics for sample "+r1_base_filename+" finished successfully"
	else:
		print "collectalignmentsummarymetrics for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd


# BUILDBAMINDEX
# =====================================================================
#BuildBamIndex I=input.bam 
bam_index= picard_output_dir+"markduplicates/"+r1_base_filename+ "_marked_duplicates.bai"
if(not os.path.isfile(bam_index)) and (("buildbamindex" in steps_to_process) or ("all" in steps_to_process)):
	cmd = BuildBamIndex[:]
	cmd.append("I=")
	cmd.append(marked_duplicates_bam)
	cmd.append("O=")
	cmd.append(bam_index)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "buildbamindex for sample "+r1_base_filename+" finished successfully"
	else:
		print "buildbamindex for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd
	subprocess.call(["mkdir", gatk_output_dir])

"""
#RealignerTargetCreator and IndelRealigner are not necessary if using Haplotypecaller for variant calling!

# REALIGNERTARGETCREATOR
# =====================================================================
#RealignerTargetCreator -T=RealignerTargetCreator -R=reference.fasta -I=input.bam --known=indels.vcf -o=forIndelRealigner.intervals
intervals_for_IR = gatk_output_dir+r1_base_filename+"_forIndelRealigner.intervals"
indels_vcf_Mills = "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
indels_vcf_1000g = "1000G_phase1.indels.hg19.sites.vcf"
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
		print "realignertargetcreator for sample "+r1_base_filename+" finished successfully"
	else:
		print "realignertargetcreator for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd
	
#RealignerTargetCreator and IndelRealigner are not necessary if using Haplotypecaller for variant calling!

# INDELREALIGNER
# =====================================================================
#IndelRealigner -T=IndelRealigner -R=reference.fasta -I=input.bam -known=indels.vcf -targetIntervals=intervalListFromRTC.intervals -o=realignedBam.bam
realigned_bam = gatk_output_dir+r1_base_filename+"_realigned.bam"
mills_intervals = "/dataVolume/storage/rb_pipeline/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
thousandg_intervals = "/dataVolume/storage/rb_pipeline/1000G_phase1.indels.hg19.sites.vcf" 
if(not os.path.isfile(realigned_bam)) and (("indelrealigner" in steps_to_process) or ("all" in steps_to_process)):
	cmd = IndelRealigner[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I")
	cmd.append(marked_duplicates_bam)
	cmd.append("-known")
	cmd.append(mills_intervals)
	cmd.append("-known")
	cmd.append(thousandg_intervals)
	cmd.append("-targetIntervals")
	cmd.append(intervals_for_IR)
	cmd.append("-o")
	cmd.append(realigned_bam)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "indelrealigner for sample "+r1_base_filename+" finished successfully"
	else:
		print "indelrealigner for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd
"""

# BASERECALIBRATOR
# =====================================================================
#BaseRecalibrator not necessary for samples to be run through mutect2 or haplotypecaller -T=IndelRealigner -R=reference.fasta -I=input.bam -known=indels.vcf -targetIntervals=intervalListFromRTC.intervals -o=realignedBam.bam 
recal_report = gatk_output_dir+r1_base_filename+"_recal_report.grp"
realigned_bam = gatk_output_dir+r1_base_filename+"_realigned.bam"
indels_vcf_dbsnp = "/dataVolume/storage/kooi_rb_pipeline/dbsnp_138.hg19.excluding_sites_after_129.vcf"
indels_vcf_Mills = "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
indels_vcf_1000g = "1000G_phase1.indels.hg19.sites.vcf"
if(not os.path.isfile(recal_report)) and (("baserecalibrator" in steps_to_process) or ("all" in steps_to_process)):
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
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "baserecalibrator for sample "+r1_base_filename+" finished successfully"
	else:
		print "baserecalibrator for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd
	
# PRINTREADS
# =====================================================================
#PrintReads -T=PrintReads -R=reference.fasta -I=input.bam -BQSR= recalibration_report.grp -o=output.bam
recal_bam = gatk_output_dir+r1_base_filename+"_recalibrated.bam"
if (not os.path.isfile(recal_bam)) and (("printreads" in steps_to_process) or ("all" in steps_to_process)):
	cmd = PrintReads[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I")
	cmd.append(realigned_bam)
	cmd.append("-BQSR")
	cmd.append(recal_report)
	cmd.append("-o")
	cmd.append(recal_bam)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "printreads for sample "+r1_base_filename+" finished successfully"
	else:
		print "printreads for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd
	
# BUILDBAMINDEX
# =====================================================================
#BuildBamIndex I=input.bam 
recal_bam_index= gatk_output_dir+r1_base_filename+"_recalibrated.bai"
if(not os.path.isfile(recal_bam_index)) and (("buildbamindex" in steps_to_process) or ("all" in steps_to_process)):
	cmd = BuildBamIndex[:]
	cmd.append("I=")
	cmd.append(recal_bam)
	cmd.append("O=")
	cmd.append(recal_bam_index)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "buildbamindex for sample "+r1_base_filename+" finished successfully"
	else:
		print "buildbamindex for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd

"""
# ANALYZECOVARIATES
# =====================================================================
#AnalyzeCovariates -T=AnalyzeCovariates -R=reference.fasta -I=input.bam -BQSR= recalibration_report.grp -o=output.bam
recal_plot = gatk_output_dir+r1_base_filename+"_recalQC.pdf"
if(not os.path.isfile(recal_plot)) and (("analyzecovariates" in steps_to_process) or ("all" in steps_to_process)):
	cmd = AnalyzeCovariates[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-before")
	cmd.append(realigned_bam)
	cmd.append("-after")
	cmd.append(recal_report)
	cmd.append("-plots")
	cmd.append(recal_plot)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "analyzecovariates for sample "+r1_base_filename+" finished successfully"
	else:
		print "analyzecovariates for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd
"""
	
# HAPLOTYPECALLER
# =====================================================================
if not os.path.isdir(haplocaller_output_dir):
	subprocess.call(["mkdir", haplocaller_output_dir])
haplocaller_gvcf = haplocaller_output_dir+r1_base_filename+"_haplotype.g.vcf"

if(not os.path.isfile(haplocaller_gvcf)) and (("haplotypecaller" in steps_to_process) or ("all" in steps_to_process)):
	cmd  = HaplotypeCaller[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-I")
	cmd.append(recal_bam)
	cmd.append("--emitRefConfidence")
	cmd.append("GVCF")
	cmd.append("--dbsnp")
	cmd.append(indels_vcf_dbsnp)
	cmd.append("-L")
	cmd.append(target_intervals)
	cmd.append("-o")
	cmd.append(haplocaller_gvcf)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "haplotypecaller for sample "+r1_base_filename+" finished successfully"
	else:
		print "haplotypecaller for sample "+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd

# GENOTYPEGVCFS
# =====================================================================
pre_filter_all_vcf = haplocaller_output_dir+"pre_filter_all.vcf"
if(not os.path.isfile(pre_filter_all_vcf)) and (("genotypegvcfs" in steps_to_process) or ("all" in steps_to_process)):
	with open('/home/thor/kooi_rb_output/gatk/genotype_gvcf_files.txt', 'r+') as f:
		for file in os.listdir(haplocaller_output_dir):
			if file.endswith(".g.vcf"):
				print(file)
		gvcf_list = f.read().splitlines()
	cmd = GenotypeGVCFs[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	for i in gvcf_list:
		cmd.append("--variant")
		cmd.append(gvcf_list[i])
	cmd.append("-stand_call_conf")
	cmd.append("50.0")
	cmd.append("-stand_emit_conf")
	cmd.append("10.0")
	cmd.append("-o")
	cmd.append(pre_filter_all_vcf)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "genotypegvcfs for sample "+r1_base_filename+" finished successfully"
	else:
		print "genotypegvcfs for sample "+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd
	

"""
# VARIANTFILTRATION
# =====================================================================
#VariantFiltration -T=VariantFiltration -R=reference.fasta -o=output.vcf --variant= input.vcf --filterExpression= "AB < 0.2 || MQ0 > 50" --filterName= "SomeFilterName"
RB_vcf = gatk_output_dir+"all_filtered_output.vcf"
if(not os.path.isfile(RB_vcf)) and (("variantfiltration" in steps_to_process) or ("all" in steps_to_process)):
	cmd = VariantFiltration[:]
	cmd.append("-R")
	cmd.append(reference_genome)
	cmd.append("-o")
	cmd.append(filtered_RB_vcf)
	cmd.append("--variant")
	cmd.append(pre_filter_all_vcf)
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "variantfiltration for sample "+r1_base_filename+" finished successfully"
	else:
		print "variantfiltration for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd


# ANNOVAR
# =====================================================================
RB_avinput = annovar_output_dir+r1_base_filename+".avinput"
if not os.path.isdir(annovar_output_dir):
	subprocess.call(["mkdir", annovar_output_dir])
#convert vcf file to annovar input file format
if(not os.path.isfile(RB_avinput)) and (("annovar" in steps_to_process) or ("all" in steps_to_process)):
	cmd = ConvertToAnnovar[:]
	cmd.append("-format")
	cmd.append("vcf4")
	cmd.append(RB_vcf)
	cmd.append("-outfile")
	cmd.append(RB_avinput)
	cmd.append("-allsample")
	cmd.append("-includeinfo")
	cmd.append("-withzyg")
	cmd.append("-withfreq")
	cmd.append("-include")
	cmd.append("-comment")
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "ConvertToAnnovar for sample "+r1_base_filename+" finished successfully"
	else:
		print "ConvertToAnnovar for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd
	
if ("annovar" in steps_to_process) or ("all" in steps_to_process):
	cmd = GeneAnnotation[:]
	cmd.append("-out")
	cmd.append(r1_base_filename)
	cmd.append("-build")
	cmd.append("hg19")
	cmd.append(RB_avinput)
	cmd.append("humandb/") #need to update?
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "GeneAnnotation for sample "+r1_base_filename+" finished successfully"
	else:
		print "GeneAnnotation for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd
	# need also region and filter-based annotations?
	#how to include information from pathogen predicting programs?

# COPYWRITER
# =====================================================================
if ("copywriter" in steps_to_process) or ("all" in steps_to_process):
	subprocess.call(["mkdir",copywriter_output_dir])
	cmd = CopyWriter[:]
	print " ".join(cmd)
	errcode = subprocess.call(cmd)
	if(errcode == 0):
		print "GeneAnnotation finished successfully"
	else:
		print "GeneAnnotation for sample"+r1_base_filename+" failed !!!!"
		del steps_to_process[:]
	del cmd

#GISTIC
# =====================================================================
if ("gistic" in steps_to_process) or ("all" in steps_to_process):
	cmd = Gistic[:]

#COMBINEVARIANTS
# =====================================================================
#combined_vcf = bafsegment_output_dir+r1_base_filename+".vcf"
#if(not os.path.isfile(combined_vcf)) and (("combinevariants" in steps_to_process) or ("all" in steps_to_process)):
#	cmd = CombineVariants[:]
#	cmd.append("-R")
#	cmd.append(reference_genome)
#	cmd.append("--variant")
#	cmd.append(vcf_inputs)

#SAMTOOLS MPILEUP
# =====================================================================
#samtools mpileup -uDl regions.bed -f hg19.fa file1.sorted.bam file2.sorted.bam | bcftools view -bvcg - > RAL_samtools.raw.bcf
#bcftools view RAL_samtools.raw.bcf | vcfutils.pl varFilter -D100 > RAL_samtools.vcf
"""
