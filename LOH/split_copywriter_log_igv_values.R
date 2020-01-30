#!/usr/bin/Rscript

library(tidyverse)

d = read.table("~/rb_pipeline/output/copywriter/vc/50kb/CNAprofiles/log2_read_counts.igv", sep="\t", skip=1, header =TRUE) 

sample_columns <- colnames(d)[grepl(".*[0-9]{2}.(CL|T)", colnames(d))]
reference_columns <- stringr::str_replace_all(sample_columns, "(CL|T)", "N")

pair_columns <- purrr::map2(sample_columns, reference_columns, ~c(.x, .y))

pair_columns <- purrr::map(pair_columns, ~purrr::set_names(.x, gsub("\\.", "-", stringr::str_extract(.x, "[0-9]{2}.(T|CL|N)"))))

names(pair_columns) <- purrr::map_chr(pair_columns, ~names(.x)[[1]])

var_columns <- c("Chromosome", "Start", "End", "Feature")

keep_list <- purrr::map(pair_columns, ~dplyr::select(d, one_of(var_columns), {{.x}}))

igv_file_paths <- fs::path("~/rb_pipeline/output/bafsegmentation/data/", paste0(names(keep_list), "_read_counts.igv"))

purrr::map2(keep_list, igv_file_paths, write_csv)

