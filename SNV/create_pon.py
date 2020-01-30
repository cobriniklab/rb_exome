#!/usr/bin/env python

import sys
import os
import subprocess
import glob
import IPython

fastq_files = sys.argv[1]
output_directory = sys.argv[2]

# with open(fastq_files) as f:
#     fastq_files = list(f)
    
with open(fastq_files) as temp_file:
  fastq_files = [line.rstrip('\n') for line in temp_file]
    
r1_base_filenames = [i.replace(".fastq.gz","").split("/")[-1] for i in fastq_files]

mutect2_output_dir=output_directory+"/mutect2/"
mutect2_pon_vcfs= [mutect2_output_dir+i+"_pon_mutect2.vcf.gz" for i in r1_base_filenames]

    
# IPython.embed()
# vcf_files = [name for name in glob.glob(vcf_dir+"*N_1_mutect2.vcf")]

mutect2_pon = mutect2_output_dir+"mutect2_pon.vcf.gz"
cmd = ["gatk", "CreateSomaticPanelOfNormals"]
for i in mutect2_pon_vcfs:
  cmd += ["-vcfs", i]
cmd += ["-O", mutect2_pon]
out_file = mutect2_output_dir +"create_pon.log"
err_file = mutect2_output_dir + "create_pon.err"
print(out_file,err_file)
print( " ".join(cmd))
with open(err_file,"w")as outerr:
    with open(out_file,"w")as outf:
        errcode = subprocess.call(cmd,stdout=outf,stderr=outerr)
        if(errcode == 0):
            print( "panel of normals created successfully")
        else:
            print( "panel of normals failed !!!!")
del cmd

