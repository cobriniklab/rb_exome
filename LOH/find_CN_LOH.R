#!/usr/bin/env Rscript


# load libraries ----------------------------------------------------------
library(tidyverse)
library(plyranges)
library(karyoploteR)


# load functions -------------------------------------------------

plot_base_seg <- function(cov_granges, name_cov_grange, tn, kp, chr_select, mid.regs, ...){
  
  tr.i <- 0.037
  tr.o <- 0.040
  kp <- kpDataBackground(kp, r0=tr.o*tn, r1=tr.o*tn+tr.i) %>%
    # kpAxis(ymin =0, y = 1, cex = 0.3, r0=(tr.o*tn), r1=(tr.o*tn+tr.i)) %>% 
    # kpPlotRegions(data = mid.regs, r0=tr.o*tn, r1=tr.o*tn+tr.i, col = NA, lty=1, lwd=0.5, border="blue", data.panel=2) %>%
    kpHeatmap(cov_granges, y = cov_granges$mBAF, ymin=0, ymax=1, r0=tr.o*tn, r1=tr.o*tn+tr.i, colors = c("blue", "white", "red")) %>% 
    kpAddLabels(labels=name_cov_grange, r0=tr.o*tn, r1=tr.o*tn+tr.i,  pos=1, label.margin = 0.035, col="red", cex=0.5) %>% 
    identity()
}

plot_all_to_file <- function(raw_cov_list, file_name, chr_select, mid.regs, ...) {
  
  png(file_name, width = 5, height = 4, units = 'in', res = 800)
  plot.params <- getDefaultPlotParams(plot.type=5)
  plot.params$ideogramheight <- 3
  plot.params$data1height <- 1
  plot.params$data1inmargin <- 1
  plot.params$data2inmargin <- 1
  plot.params$data1outmargin <- 1
  plot.params$bottommargin <- 1
  plot.params$topmargin <- 7
  plot.params$leftmargin <- 0.1
  kp <- plotKaryotype(genome="hg19", plot.type = 5, plot.params = plot.params, labels.plotter = NULL, chromosomes = chr_select) %>% 
    kpAddChromosomeNames(col="red", cex = 0.45, srt=45)
  num_baf_granges <- seq(0, length(raw_cov_list)-1)
  for (i in names(raw_cov_list)) {
    num_baf <- which(names(raw_cov_list) == i) - 1
    plot_base_seg(raw_cov_list[[i]], i, num_baf, kp, chr_select, mid.regs)
  }
  dev.off()
}

# reynolds CN-LOH ---------------------------------------------------------
AI_regions_path <- "~/rb_pipeline/doc/LOH/AI_regions_reynolds.csv"
AI_regions <- readr::read_csv(AI_regions_path) %>% 
  dplyr::mutate(Assay = gsub("-", ".", Assay)) %>% 
  identity()

AI_regions <- split(AI_regions, AI_regions$Assay) %>% 
  purrr::map(as_granges) %>% 
  identity()


reynolds_seg_path <- "~/rb_pipeline/results/SCNA/reynolds_SCNA_segments.txt"
karyo_seg <- read.table(reynolds_seg_path) %>% 
  dplyr::rename(seqnames = chrom) %>%
  # dplyr::mutate(start = as.numeric(start)) %>% 
  # dplyr::filter(seg.mean < 0) %>% 
  identity()

karyo_seg <- split(karyo_seg, karyo_seg$sample_id)

karyo_seg <- purrr::map(karyo_seg, as_granges)



test <- purrr::map2(AI_regions, karyo_seg, plyranges::setdiff_ranges)
  
test <- purrr::map2(test, AI_regions, join_overlap_inner)

# load kooi SCNA peak regions 
kooi_peak_regions <- read_csv("~/rb_pipeline/doc/SCNA/kooi_SCNA_peak_regions.csv")

kooi_peak_granges <- makeGRangesFromDataFrame(kooi_peak_regions)

baf_granges_path <- "doc/RB_exome_manuscript/LOH/CN_LOH2_reynolds.png"
plot_all_to_file(test, baf_granges_path, chr_select = "auto", kooi_peak_granges)

# VC CN-LOH ---------------------------------------------------------------
AI_regions_path <- "~/rb_pipeline/doc/LOH/AI_regions_vc.csv"
AI_regions <- readr::read_csv(AI_regions_path) %>% 
  dplyr::mutate(Assay = gsub("-", ".", Assay)) %>% 
  identity()

AI_regions <- split(AI_regions, AI_regions$Assay) %>% 
  purrr::map(as_granges) %>% 
  identity()

# output vc seg
vc_seg_path <- "~/rb_pipeline/results/SCNA/vc_SCNA_segments.txt"

karyo_seg <- read.table(vc_seg_path) %>% 
  dplyr::rename(seqnames = chrom) %>%
  # dplyr::filter(seg.mean < 0) %>%
  identity()

karyo_seg <- split(karyo_seg, karyo_seg$sample_id)

karyo_seg <- purrr::map(karyo_seg, as_granges)

test <- purrr::map2(AI_regions, karyo_seg, plyranges::setdiff_ranges)

test <- purrr::map2(test, AI_regions, join_overlap_inner)

# load kooi SCNA peak regions 
kooi_peak_regions <- read_csv("~/rb_pipeline/doc/SCNA/kooi_SCNA_peak_regions.csv")

kooi_peak_granges <- makeGRangesFromDataFrame(kooi_peak_regions)

baf_granges_path <- "doc/RB_exome_manuscript/LOH/CN_LOH_vc.png"
plot_all_to_file(test, baf_granges_path, chr_select = "auto", kooi_peak_granges)
