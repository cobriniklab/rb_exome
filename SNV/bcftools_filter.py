#!/usr/bin/env python

import sys
import os 
import subprocess

vcfs=sys.argv[1:]

for i in vcfs:
  cmd = ["bcftools", "query"]
  cmd += ["-i", r'FILTER="PASS"']
  cmd += ["-f", r'%CHROM %POS %FILTER\n']
  cmd += [i]
  out_vcf = r1_base_filename+"_mosdepth.log"
    with open(log_file,"w")as logf:
        with open(coverage_table, "w") as outf:
            errcode = subprocess.call(cmd, stdout = outf)
    del cmd


# args = ["awk", r'{OFS="\t"; print $2,$4,$5,$6}', "B3LYPD.txt"]
