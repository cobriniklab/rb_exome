#!/usr/bin/python

# import system library - takes care of reading arguments from command line
import sys
# pandas library helps with reading csv file and with statistics
import pandas as pd
import numpy as np
import IPython

#argv[1] = file with intervals and log values
#argv[2] = file with sample SNPs
#argv[3] = bin_size

sample_name = sys.argv[2].split('/')[-1]
sample_name = sample_name.replace(".txt",".")
#~ sample_name = sample_name.replace(".","")[3:]+"_"+sample_name[0:2]+"."
#~ reference_name = sample_name.split("_")[1]
#~ reference_name = "N_"+reference_name

print sample_name
#~ print reference_name 


intervals = pd.read_csv(sys.argv[1],sep=",", index_col=False, skiprows=1, low_memory=False, header=0) #names=["Chr","Start","End","Feature",sample_name+"Log R Ratio","NA"],

#~ tumor = "-T_" in intervals.columns

#~ if not tumor:
	#~ new_cols = list(intervals.columns[:-2])+[intervals.columns[-1]]+[intervals.columns[-2]]
	#~ intervals = intervals.reindex_axis(new_cols, axis=1)

intervals.columns = ["Chr","Start","End","Feature",sample_name+"Log R Ratio"]

intervals = intervals.drop(intervals.columns[[3]], axis=1)
snps = pd.read_csv(sys.argv[2],names=["Name","Chr","Position","End",sample_name+"GType",sample_name+"B Allele Freq"],sep="\t", index_col=False, skiprows=1, low_memory=False, header=0)
bin_size = int(sys.argv[3])
print bin_size
snps["Start"] = (snps["Position"].floordiv(bin_size) *bin_size)+1
snps = snps.merge(intervals, on=["Chr","Start"])
snps["Chr"] = snps["Chr"].str.replace("chr", "")

snps = snps.drop(snps.columns[[3, 4, 6, 7]], axis=1) # remove here columns you don't want in output
	
outfile = sys.argv[2].replace(".txt","_with_Rlog.csv")

snps.to_csv(outfile,sep="\t", index=False)

