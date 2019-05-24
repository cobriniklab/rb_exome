#!/usr/bin/python

import subprocess
import sys
import re
import os
import datetime
import argparse
import gzip
import pipes
import IPython

# list here all steps of the pipeline. Pipeline will run all these steps, if not requested otherwise in -s argument
steps_to_process_all = ["trimmomatic", "fastqc", "bwa", "samformatconverter", "sortsam", "removeduplicates", "buildbamindex", "baserecalibrator", "printreads", "recalbamindex", "mutect2", "mutect2_pon"]

parser = argparse.ArgumentParser(description="runs onh pipeline")
parser.add_argument("-1", "--fastq-r1", dest="f1", help="fastq R1 file REQUIRED", metavar="FILE[.gz]", required=True)
parser.add_argument("-2", "--fastq-r2", dest="f2", help="fastq R2 file REQUIRED", metavar="FILE[.gz]", required=True)
parser.add_argument("-d", "--out", dest="out_dir", help="output directory REQUIRED", metavar="DIRECTORY", required=True)
parser.add_argument("-s", "--steps-to-run", dest="steps", help="steps of pipeline to run")
parser.add_argument("-rg", "--reference", dest="reference", default= "/dataVolume/storage/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa", help="human reference genome (hg19/Grch38)", metavar="FILE[.fasta]", required=True)
parser.add_argument("--overwrite", dest="overwrite")

options = parser.parse_args()
fastq_r1_location = options.f1
fastq_r2_location = options.f2
reference_genome = options.reference
bwa_index = reference_genome.replace("WholeGenomeFasta/genome.fa", "BWAIndex/genome.fa")
output_directory = options.out_dir
output_directory_root = options.out_dir
if(options.steps == "All" or None):
    steps_to_process = steps_to_process_all
else:
    steps_to_process = options.steps.split(",")

print("Will run steps:", steps_to_process)

num_threads = "7"
# *********************************************************************
# DEFINITION OF PATHS

Trimmomatic  = ["java", "-jar", "/usr/share/java/Trimmomatic-0.36/trimmomatic-0.36.jar"]
BwaMem = ["bwa", "mem"]
SamFormatConverter = ["java",  "-jar", "/usr/share/java/picard.jar", "SamFormatConverter"]
SortSam  = ["java", "-jar", "/usr/share/java/picard.jar", "SortSam"]
MarkDuplicates = ["java", "-jar", "/usr/share/java/picard.jar", "MarkDuplicates"]
Mosdepth = ["mosdepth"]
CollectHsMetrics = ["java", "-jar", "/usr/share/java/picard.jar", "CollectHsMetrics"]
CollectIsMetrics = ["java", "-jar", "/usr/share/java/picard.jar", "CollectInsertSizeMetrics"]
CollectAlignmentSummaryMetrics = ["java", "-jar", "/usr/share/java/picard.jar", "CollectAlignmentSummaryMetrics"]
AddOrReplaceReadGroups = ["java",  "-jar", "/usr/share/java/picard.jar", "AddOrReplaceReadGroups"]
BuildBamIndex = ["java", "-jar", "/usr/share/java/picard.jar", "BuildBamIndex"]
BaseRecalibrator = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK-3.8.0.jar", "-T", "BaseRecalibrator"]
PrintReads = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK-3.8.0.jar", "-T", "PrintReads"]
FilterMutectCalls = ["gatk", "FilterMutectCalls"]
VariantAnnotator = ["java", "-jar", "/usr/share/java/GenomeAnalysisTK-3.8.0.jar", "-T", "VariantAnnotator"]
VcfAnno = ["/usr/local/bin/TOOLS/vcfanno_linux64"]
vep = ["/usr/local/bin/TOOLS/ensembl-vep/vep"]
bcftools = ["bcftools", "filter"]
Mutect2_gatk4 = ["gatk", "Mutect2"]
TableAnnovar = ["/dataVolume/storage/annovar/table_annovar.pl"]
CopyWriter = ["/media/skevin/storage/rb_pipeline/copywriter/copywriter_cobrinik.R"]


min_base_quality = "15"
min_read_length  = "25"
sliding_window   = "6:20" # window_size:min_average_quality
adapter_sequence = "AAGCAGTGGTATCAA"

fastqc_output_dir = output_directory+"fastqc/"

bwa_output_dir = output_directory+"bwa/"

picard_output_dir = output_directory+"picard/"
picard_tmp_dir = output_directory+"tmp/"
alignment_qc_output_dir = output_directory+"picard/alignment_qc/"
gatk_output_dir = output_directory + "/gatk/"
mutect2_output_dir=output_directory+"/mutect2/"
mutect_output_dir=output_directory+"/gatk/mutect/"  
annovar_output_dir = output_directory + "annovar/"
r1_base_filename = fastq_r1_location.replace(".fastq.gz","").split("/")[-1]
r2_base_filename = fastq_r2_location.replace(".fastq.gz","").split("/")[-1]
fastq_r1 = fastq_r1_location.split("/")[-1]
fastq_r2 = fastq_r2_location.split("/")[-1]
RBID = r1_base_filename.split("_")[0]
bwa_sam= r1_base_filename + ".sam"
picard_bam = output_directory+r1_base_filename + ".bam"
sortsam_bam = picard_output_dir+ r1_base_filename + "_sorted.bam"

# /DEFINITION OF PATHS
# *********************************************************************

# READ GROUP PARSING
# =====================================================================
#if("rg_parse" in steps_to_process) or ("all" in steps_to_process):
with gzip.open(fastq_r1_location) as fasta:
    s = fasta.readline().decode().split(":")
    RGPL = "Illumina"
    RGPU = "7"
    RGLB = RGPU
    s[3] = RGPU
    RGID =  "_".join(s[2:5])
    RGSM = r1_base_filename
    read_groups = "@RG\\tID:"+RGID+"\\tSM:"+RGSM+"\\tPL:"+RGPL


#~ # TRIMMOMATIC
#~ # =====================================================================
trimmomatic_output_dir = output_directory + "trimmomatic/"
fastq_r1_trimmed = trimmomatic_output_dir+ r1_base_filename + "_trimmed.fastq"
fastq_r2_trimmed = trimmomatic_output_dir+ r2_base_filename + "_trimmed.fastq"

if("trimmomatic" in steps_to_process):
    if not os.path.isdir(trimmomatic_output_dir):
        subprocess.call(["mkdir", trimmomatic_output_dir])
    print("running trimmomatic")
    cmd = Trimmomatic[:]
    #print cmd
    cmd += ["PE"]
    cmd += ["-threads", str(num_threads)]
    cmd += [fastq_r1_location, fastq_r2_location]
    cmd += [fastq_r1_trimmed, fastq_r2_trimmed]
    cmd += ["LEADING:"+min_base_quality]
    cmd += ["TRAILING:"+min_base_quality]
    cmd += ["SLIDINGWINDOW:"+sliding_window]
    cmd += ["MINLEN:"+min_read_length]
    print( " ".join(cmd))
    subprocess.call(cmd)
    del cmd

# FASTQC
# =====================================================================
fastqc_output_dir = output_directory + "fastqc/"
if not os.path.isdir(fastqc_output_dir):
    subprocess.call(["mkdir", fastqc_output_dir])
if ("fastqc" in steps_to_process):
    cmd = ["fastqc",fastq_r1_location,fastq_r2_location]
    cmd += ["-o", fastqc_output_dir]
    print( " ".join(cmd))
    errcode = subprocess.call(cmd)
    if(errcode == 0):
        print( "fastqc finished successfully")
    else:
        print( "fastqc failed !!!!")
        del steps_to_process[:]
    del cmd

# BWA ALIGNER
# =====================================================================
#BWA MEM -aM -t 6 -R '@RG' reference_genome.fa fastq_r1 fastq_r2 > file.sam 
bwa_sam= bwa_output_dir+r1_base_filename + ".sam"
if not os.path.isdir(bwa_output_dir):
    subprocess.call(["mkdir", bwa_output_dir])
if ("bwa" in steps_to_process):
    print( "beginning alignment")
    print( fastq_r1, fastq_r2)
    cmd = BwaMem[:]
    cmd.append("-M")
    cmd.append("-t")
    cmd.append("6")
    cmd.append("-R")
    cmd.append(read_groups)
    #cmd.append(log)
    cmd.append(bwa_index)
    cmd.append(fastq_r1_location) 
    cmd.append(fastq_r2_location)
    print( " ".join(cmd))
    f = open(bwa_sam, "w")
    errcode = subprocess.call(cmd, stdout=f)
    if(errcode == 0):
        print( "BWA finished successfully")
    else:
        print( "BWA failed !!!!")
        del steps_to_process[:]
    f.close()
    del cmd

# SAMFORMATCONVERTER
# =====================================================================
#SamFormatConverter I=bwa_output.sam O=converted_bwa_output.bam
picard_bam = picard_output_dir+r1_base_filename + ".bam"
if not os.path.isdir(picard_output_dir):
    subprocess.call(["mkdir", picard_output_dir])
if ("samformatconverter" in steps_to_process):
    cmd = SamFormatConverter[:]
    cmd.append("I=")
    cmd.append(bwa_sam)
    cmd.append("O=")
    cmd.append(picard_bam)
    cmd.append("MAX_RECORDS_IN_RAM=200000")
    cmd.append("TMP_DIR=")
    cmd.append(picard_tmp_dir)
    print( " ".join(cmd))
    errcode = subprocess.call(cmd)
    if(errcode == 0):
        print( "samformatconverter finished successfully")
    else:
        print( "samformatconverter failed !!!!")
        del steps_to_process[:]
    del cmd

# SORTSAM
# =====================================================================
#SortSam I=input.bam O=sorted.bam SORT_ORDER=coordinate 
sortsam_bam = picard_output_dir+r1_base_filename + "_sorted.bam"
if ("sortsam" in steps_to_process):
    print( "sortsam_bam:",sortsam_bam)
    cmd = SortSam[:]
    cmd.append("I=")
    cmd.append(picard_bam)
    cmd.append("O=")
    cmd.append(sortsam_bam)
    cmd.append("SORT_ORDER=coordinate")
    print( " ".join(cmd))
    errcode = subprocess.call(cmd)
    if(errcode == 0):
        print( "sortsam finished successfully")
    else:
        print( "sortsam failed !!!!")
        del steps_to_process[:]
        os.remove(bwa_sam)
    del cmd


# ADD OR REPLACE READ GROUPS
# =====================================================================
rg_fix_bam = picard_output_dir+r1_base_filename + "_sorted_rg.bam"

if("add_rg" in steps_to_process):
    if not os.path.isdir(picard_output_dir):
        subprocess.call(["mkdir", picard_output_dir])
    print( "running AddOrReplaceReadGroups")
    cmd = AddOrReplaceReadGroups[:]
    #print( cmd
    cmd += ["I=", sortsam_bam]
    cmd += ["O=", rg_fix_bam]
    cmd += ["RGID=", RGID]
    cmd += ["RGLB=", RGLB]
    cmd += ["RGPL=", RGPL]
    cmd += ["RGPU=", RGPU]
    cmd += ["RGSM=", RGSM] 
    print( " ".join(cmd))
    subprocess.call(cmd)
    del cmd

# MARKDUPLICATES (and remove)
# =====================================================================
#MarkDuplicates I=input.bam O=marked_duplicates.bam M=marked_dup_metrics.txt
removed_duplicates_bam= picard_output_dir + r1_base_filename + "_removed_duplicates.bam"
remove_dup_metrics = picard_output_dir+r1_base_filename+"_removed_dup_metrics.txt"
picard_tmp_dir = picard_output_dir+"tmp/"
if("removeduplicates" in steps_to_process):
    cmd = MarkDuplicates[:]
    cmd.append("I=")
    cmd.append(sortsam_bam)
    cmd.append("O=")
    cmd.append(removed_duplicates_bam)
    cmd.append("M=")
    cmd.append(remove_dup_metrics)
    cmd.append("TMP_DIR=")
    cmd.append(picard_tmp_dir)
    cmd.append("REMOVE_DUPLICATES=true")
    print( " ".join(cmd))
    errcode = subprocess.call(cmd)
    if(errcode == 0):
        print( "markduplicates finished successfully")
    else:
        print( "markduplicates failed !!!!")
        del steps_to_process[:]
    del cmd

# MOSDEPTH COVERAGE
# =====================================================================
#MOSDEPTH
coverage_table= picard_output_dir+r1_base_filename + "_coverage.txt"
coverage_dist= picard_output_dir+r1_base_filename + "_cumul.dist"
target_bed = "./bin/agilent_sureselect_v5_regions_corrected.bed"

if("mosdepth" in steps_to_process):
    cmd = Mosdepth[:]
    cmd += ["--distribution", coverage_dist]
    cmd += ["--by", target_bed]
    cmd += [removed_duplicates_bam]
    print( " ".join(cmd))
    log_file = picard_output_dir +r1_base_filename+"_mosdepth.log"
    with open(log_file,"w")as logf:
        with open(coverage_table, "w") as outf:
            errcode = subprocess.call(cmd, stdout = outf)
    del cmd

# BUILDBAMINDEX
# =====================================================================
#BuildBamIndex I=input.bam 
bam_index= picard_output_dir+r1_base_filename+ "_removed_duplicates.bam.bai"
if("buildbamindex" in steps_to_process):
    cmd = BuildBamIndex[:]
    cmd.append("I=")
    cmd.append(removed_duplicates_bam)
    cmd.append("O=")
    cmd.append(bam_index)
    print( " ".join(cmd))
    out_file = picard_output_dir + r1_base_filename+"_buildbamindex.log"
    err_file = picard_output_dir + r1_base_filename+"_buildbamindex.err"
    with open(err_file,"w")as outerr:
        with open(out_file,"w")as outf:
            errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
            if(errcode == 0):
                print( "buildbamindex finished successfully")
            else:
                print( "buildbamindex failed !!!!")
                del steps_to_process[:]
        del cmd

# BASERECALIBRATOR
# =====================================================================
#BaseRecalibrator -T=IndelRealigner -R=reference.fasta -I=input.bam -known=indels.vcf -targetIntervals=intervalListFromRTC.intervals -o=realignedBam.bam 
realigned_bam = gatk_output_dir+r1_base_filename+"_realigned.bam"
recal_report = gatk_output_dir+r1_base_filename+"_recal_report.table"
indels_vcf_dbsnp = "/home/skevin/Homo_sapiens/DBSNP/dbsnp_138.hg19.excluding_sites_after_129.vcf"
indels_vcf_Mills = "/home/skevin/Homo_sapiens/MILLS_1000G/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
indels_vcf_1000g = "/home/skevin/Homo_sapiens/1000G/1000G_phase1.indels.hg19.sites.vcf"
if("baserecalibrator" in steps_to_process):
    cmd = BaseRecalibrator[:]
    cmd.append("-R")
    cmd.append(reference_genome)
    cmd.append("-I")
    cmd.append(removed_duplicates_bam)
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
    out_file = gatk_output_dir + r1_base_filename+"_baserecalibrator.log"
    err_file = gatk_output_dir + r1_base_filename+"_baserecalibrator.err"
    print( " ".join(cmd))
    with open(err_file,"w")as outerr:
        with open(out_file,"w")as outf:
            errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
            if(errcode == 0):
                print( "baserecalibrator finished successfully")
            else:
                print( "baserecalibrator failed !!!!")
                del steps_to_process[:]
        del cmd

#~ # ANALYZECOVARIATES
#~ # =====================================================================
#~ #AnalyzeCovariates -T=AnalyzeCovariates -R=reference.fasta -I=input.bam -BQSR= recalibration_report.grp -o=output.bam
#~ recal_plot = gatk_output_dir+r1_base_filename+"_recalQC.pdf"
#~ if(not os.path.isfile(recal_plot)) and (("analyzecovariates" in steps_to_process) or ("all" in steps_to_process)):
    #~ cmd = AnalyzeCovariates[:]
    #~ cmd.append("-R")
    #~ cmd.append(reference_genome)
    #~ cmd.append("-BQSR")
    #~ cmd.append(recal_report)
    #~ cmd.append("-plots")
    #~ cmd.append(recal_plot)
    #~ print( " ".join(cmd)
    #~ errcode = subprocess.call(cmd)
    #~ if(errcode == 0):
        #~ print( "analyzecovariates finished successfully"
    #~ else:
        #~ print( "analyzecovariates failed !!!!"
        #~ del steps_to_process[:]
    #~ del cmd

    
# PRINTREADS
# =====================================================================
#PrintReads -T=PrintReads -R=reference.fasta -I=input.bam -BQSR= recalibration_report.grp -o=output.bam
recal_bam = gatk_output_dir+r1_base_filename+"_recalibrated.bam"
if ("printreads" in steps_to_process):
    cmd = PrintReads[:]
    cmd.append("-R")
    cmd.append(reference_genome)
    cmd.append("-I")
    cmd.append(removed_duplicates_bam)
    cmd.append("-BQSR")
    cmd.append(recal_report)
    cmd.append("-o")
    cmd.append(recal_bam)
    out_file = gatk_output_dir + r1_base_filename+"_printreads.log"
    err_file = gatk_output_dir + r1_base_filename+"_printreads.err"
    print( " ".join(cmd))
    with open(err_file,"w")as outerr:
        with open(out_file,"w")as outf:
            errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
            if(errcode == 0):
                print( "printreads finished successfully")
            else:
                print( "printreads failed !!!!")
                del steps_to_process[:]
    del cmd
    
# recalibrated BUILDBAMINDEX
# =====================================================================
#BuildBamIndex I=input.bam 
recal_bam_index= gatk_output_dir+r1_base_filename+ "_recalibrated.bam.bai"
if("recalbamindex" in steps_to_process):
    cmd = BuildBamIndex[:]
    cmd.append("I=")
    cmd.append(recal_bam)
    cmd.append("O=")
    cmd.append(recal_bam_index)
    print( " ".join(cmd))
    errcode = subprocess.call(cmd)
    if(errcode == 0):
        print( "buildbamindex finished successfully")
    else:
        print( "buildbamindex failed !!!!")
        del steps_to_process[:]
    del cmd

# mutect2 (gatk 4)
# =====================================================================
bait_intervals = os.path.expanduser("~/Homo_sapiens/agilent_coverage_files/agilent_baits_hg19.interval_list")
target_intervals = os.path.expanduser("~/Homo_sapiens/agilent_coverage_files/agilent_targets_hg19.interval_list")
mutect2_vcf=mutect2_output_dir+r1_base_filename+"_mutect2.vcf"
mutect2_bamout=mutect2_output_dir+r1_base_filename+"_mutect2.out.bam"
mutect2_pon_vcf=mutect2_output_dir+"m2_pon.vcf.gz"

if not os.path.isdir(mutect2_output_dir):
    subprocess.call(["mkdir", mutect2_output_dir])
if "-T" in recal_bam:
    primary_recal_bam = recal_bam
    match_normal_recal_bam = re.sub('-T_1', '-N_1', primary_recal_bam)
    normal_base_filename = r1_base_filename.replace("T", "N")
elif ("-CL" in recal_bam and not re.search('\d{3}', recal_bam)):
    primary_recal_bam = recal_bam
    match_normal_recal_bam = re.sub('-CL_1', '-N_1', primary_recal_bam)
    normal_base_filename = r1_base_filename.replace("CL", "N")
else:
    steps_to_process += ["skip_mutect2"]
if ("skip_mutect2" not in steps_to_process) and (("mutect2" in steps_to_process) or ("all" in steps_to_process)):
    cmd  = Mutect2_gatk4[:]
    cmd += ["-R", reference_genome]
    cmd += ["-I", primary_recal_bam]
    cmd += ["-I", match_normal_recal_bam]
    cmd += ["-tumor", r1_base_filename]
    cmd += ["-normal", normal_base_filename]
    cmd += ["-pon", mutect2_pon_vcf]
    cmd += ["--germline-resource", indels_vcf_dbsnp]
    # cmd += ["--max_alt_allele_in_normal_fraction", "0.10"]
    # cmd += ["--max_alt_alleles_in_normal_count", "10000000"]
    # ~ cmd += ["-bamout", mutect2_bamout]
    cmd += ["-L", target_intervals] #need to evaluate 
    cmd += ["-O", mutect2_vcf]
    out_file = mutect2_output_dir + r1_base_filename+"_mutect2.log"
    err_file = mutect2_output_dir + r1_base_filename+"_mutect2.err"
    print( " ".join(cmd))
    print(out_file,err_file)
    with open(err_file,"w")as outerr:
        with open(out_file,"w")as outf:
            errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
            if(errcode == 0):
                print( "mutect2 finished successfully")
            else:
                print( "mutect2 failed !!!!")
                del steps_to_process[:]
    del cmd

# CREATE PANEL OF NORMALS--MUTECT2 (gatk 4)
# =====================================================================
rb_intervals = "bin/rb.bed"
mutect2_pon_vcf=mutect2_output_dir+"m2_pon.vcf.gz"
mutect2_pon_bamout=mutect2_output_dir+r1_base_filename+"_mutect2.bamout"

if not os.path.isdir(mutect2_output_dir):
    subprocess.call(["mkdir", mutect2_output_dir])
if ((("-N" in recal_bam) or (re.search('\d{3}', recal_bam) is not None)) and ("mutect2_pon" in steps_to_process)):
    cmd = Mutect2_gatk4[:]
    cmd += ["-R", reference_genome]
    cmd += ["-I", recal_bam]
    cmd += ["-tumor", r1_base_filename]
    # ~ cmd += ["--germline-resource", indels_vcf_dbsnp]
    cmd += ["--disable-read-filter", "MateOnSameContigOrNoMappedMateReadFilter"]
    # ~ cmd += ["-bamout", mutect2_bamout]
    cmd += ["-L", target_intervals] #need to evaluate 
    cmd += ["-O", mutect2_pon_vcf]
    out_file = mutect2_output_dir + r1_base_filename+"_mutect2.log"
    err_file = mutect2_output_dir + r1_base_filename+"_mutect2.err"
    print(out_file,err_file)
    print( " ".join(cmd))
    with open(err_file,"w")as outerr:
        with open(out_file,"w")as outf:
            errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
            if(errcode == 0):
                print( "panel of normals finished successfully")
            else:
                print( "panel of normals failed !!!!")
                del steps_to_process[:]
    del cmd

# MANTA
# =====================================================================
bait_intervals = "./agilent_baits_hg19.interval_list"
target_intervals = "./agilent_targets_hg19.interval_list"
manta_cell_dir=output_directory+"manta/"+r1_base_filename+"/"
manta_indels=manta_cell_dir+"results/variants/candidateSmallIndels.vcf.gz"
manta_runscript=manta_cell_dir+"runWorkflow.py"
RunManta = [manta_runscript, "-m", "local", "-j", "8", "--quiet"]

if not os.path.isdir(manta_cell_dir):
    subprocess.call(["mkdir", "-p", manta_cell_dir])
if "-T" in recal_bam:
    primary_recal_bam = recal_bam
    match_normal_recal_bam = re.sub('-T_1', '-N_1', primary_recal_bam)

elif "-CL" in recal_bam:
    primary_recal_bam = recal_bam
    match_normal_recal_bam = re.sub('-CL_1', '-N_1', primary_recal_bam)

else:
    steps_to_process += "skip_manta"
if((not "skip_manta" in steps_to_process) and ("manta" in steps_to_process or "strelka" in steps_to_process)):
    cmd = ["/home/skevin/TOOLS/manta-1.2.1.centos6_x86_64/bin/configManta.py"]
    cmd += ["--normalBam="+match_normal_recal_bam]
    cmd += ["--tumorBam="+primary_recal_bam]
    cmd += ["--exome"]
    cmd += ["--referenceFasta="+reference_genome]
    cmd += ["--runDir="+manta_cell_dir] 
    print( " ".join(cmd))
    print( manta_runscript)
    p = subprocess.call(cmd)
    del cmd
    cmd = RunManta[:]
    print( " ".join(cmd))
    subprocess.call(cmd)
    del cmd


# STRELKA
# =====================================================================
strelka_output_dir=output_directory+"strelka/"
strelka_cell_dir=strelka_output_dir+r1_base_filename+"/"
manta_indels=manta_cell_dir+"results/variants/candidateSmallIndels.vcf.gz"
strelka_snv_vcf=strelka_cell_dir+"results/variants/somatic.snvs.vcf.gz" 
strelka_indel_vcf=strelka_cell_dir+"results/variants/somatic.indels.vcf.gz"

RunStrelka = [strelka_cell_dir+"runWorkflow.py", "-m", "local", "-j", "8", "--quiet"]

if not os.path.isdir(strelka_cell_dir):
    subprocess.call(["mkdir", "-p", strelka_cell_dir])
if "-T" in recal_bam:
    primary_recal_bam = recal_bam
    match_normal_recal_bam = re.sub('-T_1', '-N_1', primary_recal_bam)
elif "-CL" in recal_bam:
    primary_recal_bam = recal_bam
    match_normal_recal_bam = re.sub('-CL_1', '-N_1', primary_recal_bam)
else:
    steps_to_process += "skip_strelka"
if((not "skip_strelka" in steps_to_process) and ("strelka_call" in steps_to_process or "strelka" in steps_to_process)):
    cmd = ["/home/skevin/TOOLS/strelka-2.8.4.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py"]
    cmd += ["--normalBam", match_normal_recal_bam]
    cmd += ["--tumorBam", primary_recal_bam]
    cmd += ["--exome"]
    cmd += ["--referenceFasta", reference_genome]
    cmd += ["--indelCandidates", manta_indels]
    cmd += ["--runDir", strelka_cell_dir] 
    print( " ".join(cmd))
    subprocess.call(cmd)
    print( datetime.datetime.today())
    del cmd
    cmd = RunStrelka[:]
    print( " ".join(cmd))
    subprocess.call(cmd)
    del cmd
    
    
# FILTER MUTECT CALLS
# =====================================================================
#VariantFiltration -T=VariantFiltration -R=reference.fasta -o=output.vcf --variant= input.vcf --filterExpression= "AB < 0.2 || MQ0 > 50" --filterName= "SomeFilterName"
filtered_vcf = mutect2_output_dir+r1_base_filename+"_filtered.vcf"
if (("-N" in recal_bam) or (re.search('\d{3}', recal_bam) is not None)):
    mutect2_vcf = mutect2_pon_vcf
if("filtermutectcalls" in steps_to_process):
    cmd = FilterMutectCalls[:]
    cmd += ["-V", mutect2_vcf]
    cmd += ["-O", filtered_vcf]
    print( " ".join(cmd))
    errcode = subprocess.call(cmd)
    if(errcode == 0):
        print( "filtermutectcalls finished successfully")
    else:
        print( "filtermutectcalls failed !!!!")
        del steps_to_process[:]
    del cmd    

 # GATK ANNOTATE VARIANTS
# =====================================================================
annotated_vcf = mutect2_output_dir+r1_base_filename+"_annotated.vcf"

if("annotate_vars" in steps_to_process):
    cmd = VariantAnnotator[:]
    cmd.append("-R")
    cmd.append(reference_genome)
    cmd += ["-o", annotated_vcf]
    cmd += ["-V", filtered_vcf]
    cmd += ["-A", "VariantType"]
    cmd += ["-L", target_intervals]
    cmd += ["--dbsnp", indels_vcf_dbsnp]
    print( " ".join(cmd))
    errcode = subprocess.call(cmd)
    if(errcode == 0):
        print( "de novo annotation for sample "+r1_base_filename+" finished successfully")
    else:
        print( "de novo annotation for sample "+r1_base_filename+" failed !!!!")
    del cmd

# ~ # ANNOVAR
# ~ # =====================================================================
# ~ annovar_output_dir = output_directory+"annovar/"
# ~ annovar_vcf = annovar_output_dir+r1_base_filename+"_annovar"

# ~ if not os.path.isdir(annovar_output_dir):
    # ~ subprocess.call(["mkdir", annovar_output_dir])
# ~ if("annovar" in steps_to_process):

# ~ # table_annovar.pl example/ex1.avinput humandb/ -buildver hg19 -out myanno -remove -protocol refGene -operation g -nastring .
    # ~ cmd = TableAnnovar[:]
    # ~ cmd.append(annotated_vcf)
    # ~ cmd.append("/dataVolume/storage/annovar/humandb/")
    # ~ cmd.append("-buildver")
    # ~ cmd.append("hg19")
    # ~ cmd.append("--out")
    # ~ cmd.append(annovar_vcf)
    # ~ cmd.append("-remove")
    # ~ cmd.append("-protocol")
    # ~ cmd += ["refGene"]
    # ~ cmd.append("-operation")
    # ~ cmd += ["g"]
    # ~ cmd.append("-nastring")
    # ~ cmd.append(".")
    # ~ cmd.append("-otherinfo")
    # ~ cmd.append("-vcfinput")
    # ~ print( " ".join(cmd)
    # ~ errcode = subprocess.call(cmd)
    # ~ if(errcode == 0):
        # ~ print( "annovar finished successfully"
    # ~ else:
        # ~ print( "annovar failed !!!!"
        # ~ del steps_to_process[:]
    # ~ del cmd

      
# Run VEP on Mutect2 output
# ====================================================================
mutect2_vep_vcf = mutect2_output_dir+r1_base_filename+"_mutect2_vep.vcf"
if("mutect2_vep" in steps_to_process):
    cmd = vep[:]
    cmd += ["--cache", "--port", "3337"]
    cmd += ["--hgvs"]
    # cmd += ["--transcript_filter", '"stable_id match N[M]_"']
    cmd += ["--vcf", "--pick_allele_gene", "--exclude_predicted", "--force_overwrite"]
    cmd += ["--check_existing", "--af_gnomad"]
    cmd += ["-i", annotated_vcf]
    cmd += ["-o", mutect2_vep_vcf]
    print( " ".join(cmd))
    log_file = mutect2_output_dir +r1_base_filename+"_vcf_vep.log"
    with open(mutect2_vep_vcf,"wb") as outf, open(log_file, "wb") as outerr:
        subprocess.call(cmd, stdout=outf, stderr=outerr)
    #~ with open(log_file, "wb") as outerr:
        #~ gzed = subprocess.Popen(["gzip", ">", mutect2_vep_anno_vcf], stdin=subprocess.PIPE, stderr = outerr)
        #~ vep = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=outerr)
        #~ vep.stdout.close()
        #~ vep.wait()
        #~ gzed.wait()
        #~ gzed.communicate()
    print( datetime.datetime.today())
    del cmd
    
# MUTECT2 VCF ANNO
# =====================================================================
vcfanno_lua = os.path.expanduser("~/rb_pipeline/bin/gnomad.lua")
vcfanno_toml = os.path.expanduser("~/rb_pipeline/bin/gnomad.toml")
mutect2_anno_vcf = mutect2_output_dir+r1_base_filename+"_anno.vcf"
if("mutect2_vcf_anno" in steps_to_process):
    cmd = VcfAnno[:]
    cmd += ["-lua", vcfanno_lua]
    cmd += [vcfanno_toml]
    cmd += [mutect2_vep_vcf]
    print( " ".join(cmd))
    log_file = mutect2_output_dir +r1_base_filename+"_vcf_anno.log"
    with open(log_file,"w")as logf:
        with open(mutect2_anno_vcf, "w") as outf:
            errcode = subprocess.call(cmd, stdout = outf)
        del cmd

# STRELKA SNV VCF ANNO
# =====================================================================
strelka_snv_anno_vcf = strelka_output_dir+r1_base_filename+"_strelka_snv_anno.vcf"
if("strelka_snv_vcf_anno" in steps_to_process or "strelka" in steps_to_process):
    cmd = VcfAnno[:]
    cmd += ["-lua", vcfanno_lua]
    cmd += [vcfanno_toml]
    cmd += [strelka_snv_vcf]
    print( " ".join(cmd))
    log_file = strelka_output_dir +r1_base_filename+"_vcf_anno.log"
    with open(strelka_snv_anno_vcf,"wb") as outf, open(log_file, "wb") as outerr:
        subprocess.call(cmd, stdout=outf, stderr=outerr)
    print( datetime.datetime.today())
    del cmd
    
# STRELKA INDEL VCF ANNO
# =====================================================================
strelka_indel_anno_vcf = strelka_output_dir+r1_base_filename+"_strelka_indel_anno.vcf"
if("strelka_indel_vcf_anno" in steps_to_process or "strelka" in steps_to_process):
    cmd = VcfAnno[:]
    cmd += ["-lua", vcfanno_lua]
    cmd += [vcfanno_toml]
    cmd += [strelka_indel_vcf]
    print( " ".join(cmd))
    log_file = strelka_output_dir +r1_base_filename+"_indel_vcf_anno.log"
    with open(strelka_indel_anno_vcf,"wb") as outf, open(log_file, "wb") as outerr:
        subprocess.call(cmd, stdout=outf, stderr=outerr)
    print( datetime.datetime.today())
    del cmd
    
# Run VEP on Strelka SNV
# ====================================================================
strelka_vep_anno_vcf = strelka_output_dir+r1_base_filename+"_strelka_snv_vep.vcf"
if("strelka_snv_vep" in steps_to_process or "strelka" in steps_to_process):
    cmd = vep[:]
    cmd += ["--cache", "--port", "3337"]
    cmd += ["--hgvs"]
    # cmd += ["--transcript_filter", '"stable_id match N[M]_"']
    cmd += ["--vcf", "--pick_allele_gene", "--exclude_predicted", "--force_overwrite"]
    cmd += ["-i", strelka_snv_anno_vcf]
    cmd += ["-o", strelka_vep_anno_vcf]
    print( " ".join(cmd))
    log_file = strelka_output_dir +r1_base_filename+"_snv_vcf_vep.log"
    with open(strelka_vep_anno_vcf,"wb") as outf, open(log_file, "wb") as outerr:
        subprocess.call(cmd, stdout=outf, stderr=outerr)
    #~ with open(log_file, "wb") as outerr:
        #~ gzed = subprocess.Popen(["gzip", ">", strelka_vep_anno_vcf], stdin=subprocess.PIPE, stderr = outerr)
        #~ vep = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=outerr)
        #~ vep.stdout.close()
        #~ vep.wait()
        #~ gzed.wait()
        #~ gzed.communicate()
    print( datetime.datetime.today())
    del cmd
    
# Run VEP on Strelka indel
# ====================================================================
strelka_vep_anno_vcf = strelka_output_dir+r1_base_filename+"_strelka_indel_vep.vcf"
if("strelka_indel_vep" in steps_to_process or "strelka" in steps_to_process):
    cmd = vep[:]
    cmd += ["--cache", "--port", "3337"]
    cmd += ["--hgvs"]
    # cmd += ["--transcript_filter"]
    cmd += ["--vcf", "--pick_allele_gene", "--exclude_predicted", "--force_overwrite"]
    cmd += ["-i", strelka_indel_anno_vcf]
    cmd += ["-o", strelka_vep_anno_vcf]
    print( " ".join(cmd))
    log_file = strelka_output_dir +r1_base_filename+"_indel_vcf_vep.log"
    with open(strelka_vep_anno_vcf,"wb") as outf, open(log_file, "wb") as outerr:
        subprocess.call(cmd, stdout=outf, stderr=outerr)
    #~ with gzip.open(strelka_vep_anno_vcf,"wb") as outf, open(log_file, "wb") as outerr:
        #~ vep = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=outerr)
        #~ gzed = subprocess.Popen(["gzip"], stdout=outf, stdin=vep.stdout, stderr = outerr)
        #~ vep.stdout.close()
        #~ vep.wait()
        #~ gzed.wait()
        #~ gzed.communicate()
    print( datetime.datetime.today())
    del cmd
    
# Run bcftools
# ====================================================================
bcftools_output_dir = output_directory+"bcftools/"
bcftools_filtered_vcf = bcftools_output_dir+r1_base_filename+"_bcftools_filtered.vcf"
if("bcftools_filter" in steps_to_process):
    cmd = bcftools[:]
    cmd += ["-i", """'FILTER="PASS"'"""]
    # cmd += ["-f", """'%CHROM %POS %FILTER\n'"""]
    cmd +=[strelka_vep_anno_vcf]
    print( " ".join(cmd))
    log_file = bcftools_output_dir +r1_base_filename+"_bcftools_filtered.log"
    with open(bcftools_filtered_vcf,"wb") as outf, open(log_file, "wb") as outerr:
        subprocess.call(cmd, stdout=outf, stderr=outerr)
    print( datetime.datetime.today())
    del cmd
    
    # with open('aln.bam', 'wb',0) as output_file: 
    # cmd = Popen(["samtools", "view",'-bS','aln.sam'],
    #             stdout=(output_file))
    
    
#~ ./vep --cache --dir_cache /Software/ensembl-vep/.vep  --stats_text S39_Run3.html --refseq --hgvs --fork 4 -tab --custom /Software/ensembl-vep/.vep/score.bed.gz,score,bed,exact,0 --custom /Software/ensembl-vep/.vep/EX.bed.gz,EX,bed,exact,0 --pick_allele_gene --exclude_predicted --port 3337 -i S39_Run3.recode.vcf -o S39_Run3.txt

#~ # BAFSEGMENTATION
#~ # =====================================================================
baf_output_dir = output_directory+"bafsegmentation/"
all_samples_bed = os.path.expanduser("~/rb_pipeline/output/bafsegmentation/all_samples_merged.bed")
if("bafsegmentation" in steps_to_process):
    # get path and prefix for reference bam file
    
    sample_name = r1_base_filename.replace("_1","")
    
    ref_path=removed_duplicates_bam.replace("T", "N").replace("CL", "N")
    ref_name=sample_name.replace("T", "N").replace("CL", "N")
    
    # define output for sample and reference mpileups
    sample_mpileup=baf_output_dir+ "/data/"+sample_name+".mpileup"
    ref_mpileup="output/bafsegmentation/data/"+ref_name+".mpileup"
    
    # define output for mpileup2baf for sample and reference
    sample_txt="output/bafsegmentation/extracted/"+sample_name+".txt"
    ref_txt="output/bafsegmentation/extracted/"+ref_name+".txt"
    print(ref_path, ref_name, sample_mpileup, ref_mpileup, sample_txt, ref_txt)
    if (not os.path.isfile(sample_mpileup)):
        #~ # run samtools on each sample and matching reference
        cmd = ["samtools", "mpileup"]
        cmd += ["--min-MQ", "35"]
        cmd += ["--min-BQ", "20"]
        cmd += ["-f", reference_genome]
        cmd += ["-l", all_samples_bed]
        cmd_sample = cmd + [removed_duplicates_bam]
        cmd_sample += ["-o", sample_mpileup]
        cmd_ref = cmd + [ref_path]
        cmd_ref += ["-o", ref_mpileup]
        log_file = baf_output_dir+sample_name+"_bafsegmentation.log"
        err_file = baf_output_dir+sample_name+"_bafsegmentation.err"
        print( " ".join(cmd_sample))
        with open(log_file, "wb") as outf, open(err_file, "wb") as outerr:
            subprocess.call(cmd_sample, stdout=outf, stderr=outerr)
        with open(log_file, "wb") as outf, open(err_file, "wb") as outerr:
            subprocess.call(cmd_ref, stdout=outf, stderr=outerr)
        del cmd,cmd_sample,cmd_ref
        
    if (not os.path.isfile(sample_txt)):
    #~ # run mpileup2baf on each sample
        cmd_baf = ["perl", "bin/mpileup2baf.pl"]
        log_file = baf_output_dir+sample_name+"_bafsegmentation.log"
        err_file = baf_output_dir+sample_name+"_bafsegmentation.err"
        print( " ".join(cmd_baf))
        with open(sample_mpileup, "rb") as inf, open(sample_txt, "wb") as outf, open(err_file, "wb") as outerr:
            sample_mp2baf = subprocess.Popen(cmd_baf, stdin=inf, stdout=outf, stderr=outerr)
            sample_mp2baf.wait()
        with open(ref_mpileup, "rb") as inf, open(ref_txt, "wb") as outf, open(err_file, "wb") as outerr:
            ref_mp2baf = subprocess.Popen(cmd_baf, stdin=inf, stdout=outf, stderr=outerr)
            ref_mp2baf.wait()
        del cmd_baf
    
    #~ # add logr ratios to baf.txt file before running bafsegmentation
    #~ src/add_log_R_value_to_BAF.py output/copywriter/20kb/$sample_name/CNAprofiles/log2_read_counts.igv $sample_txt $bin_size
    #~ src/add_log_R_value_to_BAF.py output/copywriter/20kb/$sample_name/CNAprofiles/log2_read_counts.igv $reference_txt $bin_size
     
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

