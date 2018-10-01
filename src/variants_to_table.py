#!/usr/bin/python

import sys
import re
import os
import subprocess

#test GATK variantstotable

vcf_output_dir = "/home/thor/rb_output/vcf_analysis/"

VariantsToTable = ["java", "-jar", "/usr/share/java/GenomeAnalysisToolkit.jar", "-T", "VariantsToTable"]

# ANALYZECOVARIATES
# =====================================================================
#AnalyzeCovariates -T=AnalyzeCovariates -R=reference.fasta -I=input.bam -BQSR= recalibration_report.grp -o=output.bam
recal_plot = gatk_output_dir+r1_base_filename+"_recalQC.pdf"
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
