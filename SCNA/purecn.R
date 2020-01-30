#!/usr/bin/env Rscript


## ----load-packages-------------------------------------------------------
library(glue)
library(fs)
library(rprojroot)
library(tidyverse)
library(furrr)
plan(multiprocess)
proj_dir = rprojroot::find_root(criterion = has_file_pattern("*.Rproj"))
purecn_dir <- fs::path(proj_dir, "output", "PureCN")
reference_dir = "/dataVolume/storage/Homo_sapiens"

base_new_ext <- function(path, old_ext, new_ext){
  fs::path_file(gsub(old_ext, new_ext, path))
}

library(PureCN)
set.seed(1234)


## ------------------------------------------------------------------------
## ----examplegc-------------------------------------------------------------
reference.file <- fs::path(reference_dir, "UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa")
bed.file <- fs::path(reference_dir, "agilent_coverage_files/agilent_sureselect_v5_regions.bed")
mappability.file <- fs::path(reference_dir, "ENCODE", "mappability", "wgEncodeCrgMapabilityAlign100mer.bigWig")

intervals <- import(bed.file)
mappability <- import(mappability.file)

purecn_intervals_file <- fs::path(purecn_dir, "gc_file.txt")


## ---- eval = F-----------------------------------------------------------
## preprocessIntervals(intervals, reference.file, mappability=mappability, output.file = purecn_intervals_file)


## ------------------------------------------------------------------------
## ----examplecoverage-------------------------------------------------------
bam_files <- fs::path(proj_dir, "output", "gatk") %>% 
  dir_ls()

T_bam_files <- path_filter(bam_files, "*-T_1_recalibrated.bam") %>% 
  set_names(gsub("_1_recalibrated.bam", "", path_file(.)))

CL_bam_files <- path_filter(bam_files, "*-CL_1_recalibrated.bam") %>% 
  set_names(gsub("_1_recalibrated.bam", "", path_file(.)))

N_bam_files <- path_filter(bam_files, "*-N_1_recalibrated.bam") %>% 
  set_names(gsub("_1_recalibrated.bam", "", path_file(.)))


## ---- eval = F-----------------------------------------------------------
## 
## furrr::future_map(T_bam_files, ~calculateBamCoverageByInterval(.x, interval.file = purecn_intervals_file, output.file = fs::path(purecn_dir, gsub(".bam", "_coverage.txt", fs::path_file(.x)))))
## 
## furrr::future_map(CL_bam_files, ~calculateBamCoverageByInterval(.x, interval.file = purecn_intervals_file, output.file = fs::path(purecn_dir, gsub(".bam", "_coverage.txt", fs::path_file(.x)))))
## 
## furrr::future_map(N_bam_files, ~calculateBamCoverageByInterval(.x, interval.file = purecn_intervals_file, output.file = fs::path(purecn_dir, gsub(".bam", "_coverage.txt", fs::path_file(.x)))))
## 


## ------------------------------------------------------------------------
## ----example_files, message=FALSE, warning=FALSE, results='hide'-----------
tumor.coverage.files <- fs::path(purecn_dir, gsub(".bam", "_coverage.txt", fs::path_file(T_bam_files)))

cellline.coverage.files <<- fs::path(purecn_dir, gsub(".bam", "_coverage.txt", fs::path_file(CL_bam_files)))

normal.coverage.files <- fs::path(purecn_dir, gsub(".bam", "_coverage.txt", fs::path_file(N_bam_files)))

retrieve_cov_and_vcf <- function(bam_file){
  # browser()
  coverage_file <- fs::path(purecn_dir, base_new_ext(bam_file, ".bam", "_coverage.txt"))
  
  vcf_file <- fs::path(proj_dir, "output", "mutect2") %>% 
    dir_ls() %>% 
    path_filter(paste0("*", base_new_ext(bam_file, "_recalibrated.bam", "_mutect2.vcf")))
  
  return(list(coverage_file = coverage_file, vcf_file = vcf_file))
}

T_files <- purrr::map(T_bam_files, retrieve_cov_and_vcf)

CL_files <- purrr::map(CL_bam_files, retrieve_cov_and_vcf)

N_files <- purrr::map(N_bam_files, retrieve_cov_and_vcf)



## ------------------------------------------------------------------------
purecn_intervals_file <- fs::path(reference_dir, "PureCN", "gc_file.txt")


## ---- eval = F-----------------------------------------------------------
## ## ----figuregccorrect, fig.show='hide', fig.width=7, fig.height=7, warning=FALSE----
## 
## normal.coverage.files <- purrr::transpose(N_files)$coverage_file
## 
## furrr::future_map(normal.coverage.files, ~correctCoverageBias(.x, purecn_intervals_file,
##     output.file=fs::path(purecn_dir, base_new_ext(.x, "_coverage.txt", "_normal_loess.txt")), plot.bias=TRUE))
## 


## ---- eval = F-----------------------------------------------------------
## ## ----normaldb--------------------------------------------------------------
## normalDB <- createNormalDatabase(unlist(normal.coverage.files))
## 
## # serialize, so that we need to do this only once for each assay
## normalDB_path <- fs::path(purecn_dir, "normalDB.rds")
## saveRDS(normalDB, file=normalDB_path)


## ------------------------------------------------------------------------
## ----normaldbpca-----------------------------------------------------------
normalDB_path <- fs::path(proj_dir, "output", "PureCN", "normalDB.rds")
normalDB <- readRDS(normalDB_path)


## ---- eval = F-----------------------------------------------------------
## tumor.coverage.files <- purrr::transpose(T_files)$coverage_file
## tumor.pools <- furrr::future_map(tumor.coverage.files, calculateTangentNormal, normalDB)
## tumor.pools_file <- fs::path(purecn_dir, "tumor.pools.rds")
## saveRDS(tumor.pools, tumor.pools_file)


## ---- eval = F-----------------------------------------------------------
## cellline.coverage.files <- purrr::transpose(T_files)$coverage_file
## cellline.pools <- furrr::future_map(cellline.coverage.files, calculateTangentNormal, normalDB)
## cellline.pools_file <- fs::path(purecn_dir, "cellline.pools.rds")
## saveRDS(cellline.pools, cellline.pools_file)


## ------------------------------------------------------------------------
tumor.pools_file <- fs::path(purecn_dir, "tumor.pools.rds")
tumor.pools <- readRDS(tumor.pools_file)

cellline.pools_file <- fs::path(purecn_dir, "cellline.pools.rds")
cellline.pools <- readRDS(cellline.pools_file)


## ------------------------------------------------------------------------
## ----calculatemappingbias--------------------------------------------------
# speed-up future runtimes by pre-calculating variant mapping biases

normal_panel_vcf_file <- fs::path(proj_dir, "output", "mutect2") %>% 
  dir_ls() %>% 
  path_filter("*m2_pon.vcf.gz")

bias <- calculateMappingBiasVcf(normal_panel_vcf_file, genome = "h19")

mapping_bias_file <- fs::path(purecn_dir, "mapping_bias.rds")
saveRDS(bias, mapping_bias_file)


## ------------------------------------------------------------------------
mapping_bias_file <- fs::path(purecn_dir, "mapping_bias.rds")
bias <- readRDS(mapping_bias_file)


## ------------------------------------------------------------------------
## ----targetweightfile1, message=FALSE--------------------------------------
interval.weight.file <- fs::path(purecn_dir, "interval_weights.txt")
calculateIntervalWeights(normalDB$normal.coverage.files, interval.weight.file)


## ---- eval = F-----------------------------------------------------------
## ## ----ucsc_segmental--------------------------------------------------------
## # Instead of using a pool of normals to find low quality regions,
## # we use suitable BED files, for example from the UCSC genome browser.
## 
## # We do not download these in this vignette to avoid build failures
## # due to internet connectivity problems.
## downloadFromUCSC <- FALSE
## if (downloadFromUCSC) {
##     library(rtracklayer)
##     mySession <- browserSession("UCSC")
##     genome(mySession) <- "hg19"
##     simpleRepeats <- track( ucscTableQuery(mySession,
##         track="Simple Repeats", table="simpleRepeat"))
##     export(simpleRepeats, "hg19_simpleRepeats.bed")
## }
## 
## snp.blacklist <- "hg19_simpleRepeats.bed"


## ------------------------------------------------------------------------
run_absolute_cn <-  function(sample_id, sample.pools, sample.coverage.files, sample.vcf.files, interval.file, normalDB, interval.weight.file) {
  browser()
  
  pool =  sample.pools[[sample_id]]
  sample.coverage.file <- tumor.coverage.files[[sample_id]]
  sample.vcf.file <- sample.vcf.files[[sample_id]]
  
  runAbsoluteCN(
      normal.coverage.file = pool,
      tumor.coverage.file = sample.coverage.file,
      vcf.file = sample.vcf.file,
      genome = "hg19",
      sampleid = sample_id,
      interval.file = interval.file,
      normalDB = normalDB,
      # args.setMappingBiasVcf=list(normal.panel.vcf.file=normal.panel.vcf.file),
      # args.filterVcf=list(snp.blacklist=snp.blacklist,
      # stats.file=mutect.stats.file),
      args.segmentation = list(interval.weight.file = interval.weight.file),
      post.optimize = FALSE,
      plot.cnv = FALSE,
      verbose = FALSE
    )
}


## ------------------------------------------------------------------------
tumor_ret <- furrr::future_map(names(T_bam_files), run_absolute_cn, 
           sample.pools = tumor.pools, 
           sample.coverage.files = purrr::transpose(T_files)$coverage_file, 
           sample.vcf.files = purrr::transpose(T_files)$vcf_file,
           interval.file = purecn_intervals_file,
           normalDB = normalDB,
           interval.weight.file = interval.weight.file)

purecn_files <- fs::path(purecn_dir, paste0(names(tumor_ret), "_PureCN.rds"))
furrr::future_map2(tumor_ret, purecn_files, ~saveRDS(.x, filename = .y))


## ------------------------------------------------------------------------
cellline_ret <- furrr::future_map(names(CL_bam_files), run_absolute_cn, 
           sample.pools = cellline.pools, 
           sample.coverage.files = purrr::transpose(CL_files)$coverage_file, 
           sample.vcf.files = purrr::transpose(CL_files)$vcf_file,
           interval.file = purecn_intervals_file,
           normalDB = normalDB,
           interval.weight.file = interval.weight.file)

purecn_files <- fs::path(purecn_dir, paste0(names(cellline_ret), "_PureCN.rds"))
furrr::future_map2(cellline_ret, purecn_files, ~saveRDS(.x, filename = .y))

