#!/usr/bin/env Rscript

library(tidyverse)
library(stchlkExome)
suppressMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86
library(S4Vectors)
library(ensemblVEP)
library(plyranges)
library(glue)
library(data.table)
library(rprojroot)
library(fs)
proj_dir = proj_dir = rprojroot::find_root(criterion = has_file_pattern("*.Rproj"))

sample_types <- c("tumors", "cell_lines", "normals", "reynolds")


## @knitr load-vcf-lists

m2_cell_line_list <- "~/rb_pipeline/results/SNV/m2_cell_line_list.rds"
m2_cell_line_list <- readRDS(m2_cell_line_list)

m2_tumor_list <- "~/rb_pipeline/results/SNV/m2_tumor_list.rds"
m2_tumor_list <- readRDS(m2_tumor_list)

m2_pon_list <- "~/rb_pipeline/results/SNV/m2_normal_list.rds"
m2_pon_list <- readRDS(m2_pon_list)

m2_reynolds_path <- "~/rb_pipeline/results/SNV/m2_reynolds_list.rds"
m2_reynolds_list <- readRDS(m2_reynolds_path)


## @knitr collate-vcfs

m2_tidy_tumors0 <- collate_vcfs(m2_tumor_list, mutect2_tn_tidy, "tn") 

m2_tidy_cell_lines0 <- collate_vcfs(m2_cell_line_list, mutect2_tn_tidy, "tn")

m2_tidy_normals0 <- collate_vcfs(m2_pon_list, mutect2_pon_tidy, "pon") 

m2_tidy_reynolds0 <- collate_vcfs(m2_reynolds_list, mutect2_pon_tidy, "pon") 

sample_types <- c("tumors", "cell_lines", "normals", "reynolds")

# save raw variant caller files------------------------------
m2_collated_vars_ps <- fs::path(proj_dir, "results", "SNV", paste0("tidy_m2_", sample_types, ".rds"))

m2_collated_vars <- list("tumor" = m2_tidy_tumors0, "cell_line" = m2_tidy_cell_lines0, "normal" = m2_tidy_normals0, "reynolds" = m2_tidy_reynolds0)

map2(m2_collated_vars, m2_collated_vars_ps, saveRDS)





