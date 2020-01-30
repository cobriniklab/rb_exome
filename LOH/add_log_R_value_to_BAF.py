#!/usr/bin/env python

# import system library - takes care of reading arguments from command line
import sys
import os
# pandas library helps with reading csv file and with statistics
import pandas as pd
import numpy as np
import IPython
#argv[1] = file with intervals and log values
#argv[2] = file with sample SNPs
#argv[3] = file with reference SNPs
#argv[4] = bin_size

sample_name = sys.argv[2].split('/')[-1]
sample_name = sample_name.replace(".txt",".")
sample_name = sample_name.replace(".","")[3:]+"_"+sample_name[0:2]+"."
reference_name = sample_name.split("_")[1]
reference_name = "N_"+reference_name

print(sample_name)
print(reference_name)
intervals = pd.read_csv(sys.argv[1],sep=",", index_col=False, low_memory=False) #names=["Chr","Start","End","Feature",sample_name+"Log R Ratio","NA"],
tumor = "-T_" in intervals.columns
# IPython.embed()
# if not tumor:
#	new_cols = list(intervals.columns[:-2])+[intervals.columns[-1]]+[intervals.columns[-2]]
#	intervals = intervals.reindex_axis(new_cols, axis=1)

intervals.columns = ["Chr","Start","End","Feature",sample_name+"Log R Ratio", reference_name+"Log R Ratio"]


sample_snps = pd.read_csv(sys.argv[2],names=["snpid", "Chr","Position","End",sample_name+"GType",sample_name+"B Allele Freq"],sep="\t", index_col=False, skiprows=1, low_memory=False, header=0)
reference_snps = pd.read_csv(sys.argv[3],names=["snpid", "Chr","Position","End",reference_name+"GType",reference_name+"B Allele Freq"],sep="\t", index_col=False, skiprows=1, low_memory=False, header=0)
sample_snps.drop(columns=["snpid", "End"], inplace=True)
reference_snps.drop(columns=["snpid", "End"], inplace=True)
snps = sample_snps.merge(reference_snps, on=["Chr", "Position"])
# IPython.embed()
snps.loc[snps[sample_name+"B Allele Freq"] < 0.25, sample_name+"GType"] = 'AA'
snps.loc[(snps[sample_name+"B Allele Freq"] >= 0.25) & (snps[sample_name+"B Allele Freq"] < 0.75), sample_name+"GType"] = 'AB'
snps.loc[snps[sample_name+"B Allele Freq"] >= 0.75, sample_name+"GType"] = 'AA'

snps.loc[snps[reference_name+"B Allele Freq"] < 0.25, reference_name+"GType"] = 'AA'
snps.loc[(snps[reference_name+"B Allele Freq"] >= 0.25) & (snps[reference_name+"B Allele Freq"] < 0.75), reference_name+"GType"] = 'AB'
snps.loc[snps[reference_name+"B Allele Freq"] >= 0.75, reference_name+"GType"] = 'AA'

bin_size = int(sys.argv[4])
print(bin_size)
snps["Start"] = (snps["Position"].floordiv(bin_size) * bin_size)+1
# intervals["Chr"] = intervals["Chr"].str.replace("chr", "")

snps = snps.merge(intervals, on=["Chr","Start"])
snps.insert(0, 'Name', range(1, 1+ len(snps)))
snps = snps[['Name','Chr', 'Position', sample_name+'GType', sample_name+'B Allele Freq', sample_name+'Log R Ratio', reference_name+'GType', reference_name+'B Allele Freq', reference_name+'Log R Ratio']]
snps.Chr = snps.Chr.str.replace("chr","")
#snps = snps.drop(snps.columns[[7, 8, 9]], axis=1) # remove here columns you don't want in output
#new_cols = list(snps.columns[0:5])+[snps.columns[-2]]+list(snps.columns[5:7])+[snps.columns[-1]]
#snps = snps.reindex_axis(new_cols, axis=1)
outfile = sys.argv[2].replace(".txt","_with_Rlog.csv")

snps.to_csv(outfile,sep="\t", index=False)

# write extracted/normal_sample_names for paired bafsegmentation workflow
sample_extract=sample_name.replace(".","")+"_extracted.txt"
reference_extract=reference_name.replace(".","")+"_extracted.txt"
normal_sample_names_txt=os.path.expanduser("~/rb_pipeline/output/bafsegmentation/extracted/normal_sample_names.txt")

# ~ paired_names = pd.DataFrame({'FilenameAssay' : sample_extract, 'FilenameNormal' : reference_extract}, index=[0])

with open(normal_sample_names_txt, "a") as outf:
	outf.write(sample_extract+"\t"+reference_extract+"\n")
