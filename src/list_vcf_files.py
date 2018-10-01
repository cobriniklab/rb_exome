#!/usr/bin/python

import os
import sys
import re
import subprocess

vcfmerge = ["vcf-merge"]

vcf_list = []
vcf_list_path = "/dataVolume/storage/rb_pipeline/vcf_list.txt"

with open(vcf_list_path, "r") as f:
	for i in f:
		vcf_list+=i 
print vcf_list
cmd = vcfmerge[:]
for i in vcf_list:
	cmd.append(i)
print " ".join(cmd)
errcode = subprocess.call(cmd)
if(errcode == 0):
	print "markduplicates finished successfully"
else:
	print "markduplicates failed !!!!"
        
del cmd

