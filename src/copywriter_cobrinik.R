#!/usr/bin/Rscript

# load required packages
library(CopywriteR)
library(CopyhelpeR)
library(DNAcopy)
library(BiocParallel)

humanreadable <- function(bin_size){
  bin_size <- bin_size/1000
  bin_size <- paste0(toString(bin_size), "kb")
}

run_copywriter <- function(bin_size, in_dir, out_dir, baits_file){
  # browser()
  dir.create(out_dir)
  
  #preCopywriteR
  
  data.folder <- tools::file_path_as_absolute(file.path(in_dir))
  preCopywriteR(output.folder = file.path(out_dir), bin.size = bin_size, ref.genome = "hg19", "chr")
  
  #BiocParallel
  bp.param <- MulticoreParam(workers = 7)
  bp.param
  
  #CopywriteR
  path <- tools::file_path_as_absolute(file.path(in_dir))
  cell_line_pattern = "^\\d{3}-CL.*_sorted.bam$"
  tumor_pattern = "^\\d{2}-T.*_sorted.bam$"
  samples_tumors <- list.files(path = path, pattern = tumor_pattern, full.names = TRUE)
  samples_cell_lines <- list.files(path = path, pattern = cell_line_pattern, full.names = TRUE)
  tumor_controls <- gsub("T", "N", samples_tumors)
  cell_line_controls <- gsub("CL", "N", samples_cell_lines)
  
  cell_line_tst <- grep("33", samples_cell_lines, value = TRUE)
  tumor_tst <- grep("43", samples_tumors, value = TRUE)
  control_tst <- grep("43|33", cell_line_controls, value = TRUE)
  
  # tst_control = grep("33|", tumor_controls, value = TRUE)
  # tst_sample = grep("31", samples_tumors, value = TRUE)

  # tst.sample.control <- data.frame(sample = c(cell_line_tst, tumor_tst, control_tst), 
  #                              control = c(control_tst, control_tst))
    
  sample.control <- data.frame(sample = c(samples_cell_lines), 
                               control = c(samples_cell_lines))
  CopywriteR(
    sample.control = sample.control,
    destination.folder = file.path(out_dir),
    reference.folder = file.path(out_dir, paste0("hg19_", humanreadable(bin_size), "_chr")),
    bp.param = bp.param,
    capture.regions.file = baits_file,
    keep.intermediary.files = F
  )
  
  #segment and visualize results
  plotCNA(destination.folder = file.path(out_dir))
  
}


bin_size = 20000
copywriter_input_dir = "~/rb_pipeline/output/picard"
copywriter_output_dir = paste0("~/rb_pipeline/output/copywriter/", humanreadable(bin_size))
copywriter_capture_regions = "~/rb_pipeline/corrected_agilent_regions.bed" #must be a .bed file!

run_copywriter(bin_size, copywriter_input_dir, copywriter_output_dir, copywriter_capture_regions)



