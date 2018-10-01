#!/usr/bin/Rscript



# load required libraries -------------------------------------------------
library(tidyverse)
library(karyoploteR)


# load required functions -------------------------------------------------

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

plot_chrom_to_file <- function(raw_cov_list, file_name, chr_select, ...) {
  png("Plot3.png", width = 4, height = 4, units = 'in', res = 800)
  plot.params <- getDefaultPlotParams(plot.type=1)
  plot.params <- getDefaultPlotParams(plot.type=1)
  plot.params$ideogramheight <- 5
  plot.params$data1height <- 150
  plot.params$data1inmargin <- 1
  plot.params$bottommargin <- 20
  plot.params$topmargin <- 20
  kp <- plotKaryotype(genome="hg19", plot.type = 1, plot.params = plot.params, labels.plotter = NULL, chromosomes = chr_select) %>% 
    kpAddChromosomeNames(col="red", cex = 0.3) %>% 
    kpAddBaseNumbers()
  
  map2(raw_cov_list, num_baf_granges, plot_base_seg, kp, chr_select)
  dev.off()
}


plot_path = "~/rb_pipeline/output/bafsegmentation/r_processed_plots/"

plot_bafs_stacked <- function(raw_cov_data, plot_title, ...){
  file_dest <- paste0(plot_path, plot_title, ".png")
  plot.params <- getDefaultPlotParams(plot.type=1)
  plot.params$ideogramheight <- 5
  plot.params$data1height <- 10
  plot.params$data1inmargin <- 1
  plot.params$bottommargin <- 20
  plot.params$topmargin <- 20
  kp <- plotKaryotype(genome="hg19", plot.type = 1, main = plot_title, plot.params = plot.params, labels.plotter = NULL, chromosomes = "auto") %>% 
    kpDataBackground(r0=0.68, r1=0.97) %>% 
    kpAddChromosomeNames(col="red", srt=30) %>% 
    kpDataBackground() %>%  
    kpAxis(ymin =0, y = 1, cex = 1.0) %>% 
    kpHeatmap(raw_cov_data, y = raw_cov_data$mBAF, ymin=0, ymax=1) 
  
}


plot_bafs_contiguous <- function(raw_cov_data, plot_title, ...){
  file_dest <- paste0(plot_path, plot_title, ".png")
  plot.params <- getDefaultPlotParams(plot.type=4)
  plot.params$ideogramheight <- 5
  plot.params$data1height <- 10
  plot.params$data1inmargin <- 1
  plot.params$bottommargin <- 20
  plot.params$topmargin <- 20
  kp <- plotKaryotype(genome="hg19", plot.type = 4, main = plot_title, plot.params = plot.params, labels.plotter = NULL, chromosomes = "auto") %>% 
    kpDataBackground(r0=0.68, r1=0.97) %>% 
    kpAddChromosomeNames(col="red", srt=30) %>% 
    kpDataBackground() %>%  
    kpAxis(ymin =0, y = 1, cex = 1.0) %>% 
    kpHeatmap(raw_cov_data, y = raw_cov_data$mBAF, ymin=0, ymax=1) 
  
}

seg_span_gene <- function(baf_segment, chr, gene_start, gene_end){
  baf_segment <- dplyr::filter(baf_segment, Chr == chr & (!(Start < gene_start & End < gene_end) & !(Start > gene_end & End > gene_end)))
}

sort_t_b4_cl <- function(baf_granges){

  sample_ids <- names(baf_granges)
  
  sample_ids2 <- strsplit(sample_ids, "_")
  
  sample_ids2 <- lapply(sample_ids2, function(x) rev(x))
  
  sample_ids2 <- sort(unlist(lapply(sample_ids2, "[", 1)))
  
  suffixes <- rep(c("T", "CL"), length(sample_ids2)/2)
  
  match_names_order <- unlist(purrr::map2(sample_ids2, suffixes, function(x,y) paste0(y, "_", x)))
  
  new_names_order <- unlist(purrr::map2(sample_ids2, suffixes, function(x,y) paste0(x, "_", y)))
  
  baf_granges <- baf_granges[order(match(sample_ids, match_names_order))]
  names(baf_granges) <- new_names_order
  return(baf_granges)
}


#chromosome names need to be in format "chrX"
clean_baf_segment <- function(baf_segment){
  baf_segment <- mutate(baf_segment, Chr = paste0("chr", Chr)) %>% 
    mutate(Chr = ifelse(Chr == "chr23", "chrX", ifelse(Chr == "chr24", "chrY", Chr))) %>% 
    dplyr::filter(mBAF > 0.7 & mBAF <= 1.00) %>%
    dplyr::filter(!grepl("N", Assay))
}


# end functions -----------------------------------------------------------

# load kooi SCNA peak regions 
kooi_peak_regions <- read.table("~/rb_pipeline/doc/SCNA/kooi_SCNA_peak_regions.csv", stringsAsFactors = F)

kooi_peak_granges <- makeGRangesFromDataFrame(kooi_peak_regions)


# load input files when processed as group ----------------------------------------------------
# baf_segment <- read.table("~/rb_pipeline/output/bafsegmentation/segmented/AI_regions.txt", header = TRUE)
# baf_segment_list <- split(baf_segment, f = baf_segment$Assay)

if (cohort == "reynolds"){
  baf_segment_files <- list.files("~/rb_pipeline/output/bafsegmentation/segmented/", recursive = TRUE, pattern = "AI_regions.txt", full.names = TRUE)
  
  # filter for reynolds baf segment files
  baf_segment_files <- baf_segment_files[grepl("\\d{3}-CL/", baf_segment_files)]
  
  baf_segment_names <- basename(dirname(baf_segment_files))
  
  baf_segment_list <- lapply(baf_segment_files, read.table, header = TRUE)
  baf_segment_list <- lapply(baf_segment_list, clean_baf_segment)
  
  baf_granges <- map(baf_segment_list, makeGRangesFromDataFrame, keep.extra.columns = TRUE)
  
  names(baf_granges) <- sapply(baf_granges, function(x) unique(x$Assay))
  baf_granges_path <- "~/rb_pipeline/doc/LOH/LOH_reynolds_plot_coverage_ai_4.png"
  plot_all_to_file(baf_granges, baf_granges_path, chr_select = "auto", kooi_peak_granges)
  
  # output table of Allelic Imbalance Regions
  AI_regions <- rbindlist(lapply(baf_granges, data.frame))
  
  AI_regions_path <- "~/rb_pipeline/doc/LOH/AI_regions_reynolds.csv"
  write.csv(AI_regions, AI_regions_path)
  
} else if (cohort == "vc"){
  # load input files when processed per sample -------------------------------
  baf_segment_files <- list.files("~/rb_pipeline/output/bafsegmentation/segmented", recursive = TRUE, pattern = "all_AI_regions.txt", full.names = TRUE)
  
  baf_segment_files <- c(
    "/home/skevin/rb_pipeline/output/bafsegmentation/segmented/33-CL/AI_regions.txt",
    "/home/skevin/rb_pipeline/output/bafsegmentation/segmented/43-T/AI_regions.txt",
    "/home/skevin/rb_pipeline/output/bafsegmentation/segmented/all_AI_regions.txt"
  )  
  # baf_segment_files <- head(baf_segment_files, -2)
  
  baf_segment_names <- basename(dirname(baf_segment_files))
  
  # for combining separate rounds of bafsegmentation
  addl_samples <- strsplit(baf_segment_names[1:2], "-")
  addl_samples <- sapply(addl_samples, function(x) paste0(x[2], "_", x[1]))
    
  baf_segment_list <- lapply(baf_segment_files, read.table, header = TRUE)
  baf_segment_list <- lapply(baf_segment_list, clean_baf_segment)  
  
  multi_baf <- split(baf_segment_list[[3]], baf_segment_list[[3]]$Assay)
  baf_segment_list <- c(baf_segment_list[1], baf_segment_list[2], multi_baf)
  baf_segment_names[1:2] <- addl_samples
  
  baf_granges <- map(baf_segment_list, makeGRangesFromDataFrame, keep.extra.columns = TRUE)
  names(baf_granges) <- sapply(baf_granges, function(x) unique(x$Assay))
  baf_granges <- sort_t_b4_cl(baf_granges)
  
  baf_granges_path <- "~/rb_pipeline/doc/LOH/LOH_vc_plot_coverage.png"
  plot_all_to_file(baf_granges, baf_granges_path, chr_select = "auto", kooi_peak_granges)
  
  # output table of Allelic Imbalance Regions
  AI_regions <- rbindlist(lapply(baf_granges, data.frame))
  
  AI_regions_path <- "~/rb_pipeline/doc/LOH/AI_regions_vc.csv"
  write.csv(AI_regions, AI_regions_path)
  
}



# Allelic Imbalance called LOH when mBAF > ?; threshold is set at 0.7 in Kooi paper

rb1_start = 48877883
rb1_end = 49056026

mycn_start <- 16080683
mycn_end <- 16087129

rb_segments <- lapply(baf_segment_list, seg_span_gene, "chr13", rb1_start, rb1_end)

mycn_segments <- baf_rb <- lapply(baf_segment_list, seg_span_gene, "chr2", mycn_start, mycn_end)

# plot contiuous, multiple samples--heatmaps ------------------------------

plot_path = "./results/bafsegmentation/karyoploter"

check_na <- function(granges){
  table(granges$mBAF)
}

lapply(baf_granges, check_na)

names(baf_granges)



