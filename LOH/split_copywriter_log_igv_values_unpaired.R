#!/usr/bin/Rscript

library(dplyr)
library(tidyr)
library(purrr)


d = read.table("~/rb_pipeline/output/copywriter/20kb/CNAprofiles/log2_read_counts.igv", sep="\t", skip=1, header =TRUE) 

dr <- dplyr::select(d, -grep(".*1$|N", colnames(d)))

# create a list of dataframes, one per sample containing chr start end feature and sample-name-logRratio
keep_list <- lapply(colnames(dr)[5:length(colnames(dr))], function(x) dplyr::select(d, c(1:4, x)))

new_names <- sapply(keep_list, function(x) colnames(x[5]))
new_names <- gsub("\\.", "-", new_names)

fixed_names <- paste0("~/rb_pipeline/output/bafsegmentation/data/", gsub("_.*", "", new_names), "_read_counts.igv")

names(keep_list) <- sapply(keep_list, function(x) colnames(x[5]))

purrr::map2(keep_list, fixed_names, write.csv, row.names=FALSE)

