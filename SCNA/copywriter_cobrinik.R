#!/usr/bin/env Rscript

## ------------------------------------------------------------------------

# load required packages
library(CopywriteR)
library(CopyhelpeR)
library(DNAcopy)
library(BiocParallel)
library('tidyverse')
library('fs')
library('rprojroot')
library('glue')


#BiocParallel
bp.param <- MulticoreParam(workers = 7)

proj_dir = rprojroot::find_root(criterion = has_file_pattern("*.Rproj"))



## ------------------------------------------------------------------------

humanreadable <- function(bin_size){
  bin_size <- bin_size/1000
  bin_size <- paste0(toString(bin_size), "kb")
}



run_copywriter <- function(bin_size, sample_files, control_files, out_dir, baits_file, ...){
  # 
  dir.create(out_dir)
  
  #preCopywriteR
  
  # data.folder <- tools::file_path_as_absolute(file.path(in_dir))
  preCopywriteR(output.folder = file.path(out_dir), bin.size = bin_size, ref.genome = "hg19", "chr")
  

  
  if(!all(sample_files == control_files)){
    sample_files <- c(sample_files, control_files)
    
    control_files <- c(control_files, control_files)
  }
  

  sample.control <- data.frame(sample = c(sample_files), 
                               control = c(control_files))
  
  CopywriteR(
    sample.control = sample.control,
    destination.folder = file.path(out_dir),
    reference.folder = file.path(out_dir, paste0("hg19_", humanreadable(bin_size), "_chr")),
    bp.param = bp.param,
    capture.regions.file = baits_file,
    ...
  )
  
  #segment and visualize results
  plotCNA(destination.folder = file.path(out_dir))
  
}



## ---- eval = F-----------------------------------------------------------
##   #CopywriteR
## 
##   path = fs::path(proj_dir, "output", "gatk")
## 
##   cell_line_pattern = "^\\d{2}-CL.*_recalibrated.bam$"
## 
##   tumor_pattern = "^\\d{2}-T.*_recalibrated.bam$"
## 
##   exp_pattern = paste(c(cell_line_pattern, tumor_pattern), collapse = "|")
## 
##   all_normal_pattern <- "^\\d{2}-N.*_recalibrated.bam$"
## 
## 
##   tumor_samples <- list.files(path = path, pattern = tumor_pattern, full.names = TRUE)
## 
##   cell_line_samples <- list.files(path = path, pattern = cell_line_pattern, full.names = TRUE)
## 
##   exp_samples <- list.files(path = path, pattern = all_normal_pattern, full.names = TRUE)
## 
##   tumor_normals <- gsub("T", "N", tumor_samples)
##   cell_line_normals <- gsub("CL", "N", cell_line_samples)
##   exp_normals <- gsub("CL|T", "N", exp_samples)
## 
##   sample_files <- list.files(path = path, pattern = exp_pattern, full.names = TRUE)
##   control_files <- gsub("CL|T", "N", sample_files)
## 


## ---- eval = F-----------------------------------------------------------
## path = "~/rb_pipeline/output/gatk"
## 
## cl_31_pattern = "^31-CL.*_recalibrated.bam$"
## t_31_pattern = "^31-T.*_recalibrated.bam$"
## exp_pattern = paste(c(cl_31_pattern, t_31_pattern), collapse = "|")
## 
## sample_files <- list.files(path = path, pattern = exp_pattern, full.names = TRUE)
## control_files <- gsub("CL|T", "N", sample_files)
## 


## ------------------------------------------------------------------------
  
  path = fs::path(proj_dir, "output", "gatk")

  cell_line_pattern = "^\\d{3}-CL.*_recalibrated.bam$"
  
  cell_line_samples <- list.files(path = path, pattern = cell_line_pattern, full.names = TRUE)
  
  sample_files <- cell_line_samples
  control_files <- cell_line_samples


## ------------------------------------------------------------------------

bin_size = 50000
copywriter_output_dir = paste0("~/rb_pipeline/output/copywriter/", humanreadable(bin_size))
copywriter_capture_regions = "~/rb_pipeline/corrected_agilent_regions.bed" #must be a .bed file!

run_copywriter(bin_size, sample_files, control_files, copywriter_output_dir, copywriter_capture_regions, keep.intermediary.files = T)


