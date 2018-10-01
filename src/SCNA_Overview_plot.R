#!/usr/bin/Rscript

# cohort = "reynolds"

# load required libraries -------------------------------------------------


library(reshape2)
library(tidyverse)
library(data.table)
library(karyoploteR)
library(purrr)
library(scales)
library(biobroom)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene #shorthand (for convenience)
library(org.Hs.eg.db)
library(plyranges)

# load required functions --------------------------------------------------

sort_T_b4_CL <- function(karyo_granges){

  karyo_names <- split(names(karyo_granges), gsub("\\..*", "", names(karyo_granges)))
  karyo_names <- lapply(karyo_names, sort, decreasing = TRUE)
  
  karyo_granges <- karyo_granges[unlist(karyo_names)]  
}

load_seg_objs <- function(segmentation_files){
  if (length(segmentation_files) == 1){
    segmentation_objects <- get(load(segmentation_files))
    segmentation_objects <- split(segmentation_objects$output, segmentation_objects$output$ID)
    return(segmentation_objects)
  } else{   # if loading rdata files containing single sample per file ---------------
    segment_names <- strsplit(segmentation_files,'/')
    segment_names <- sapply(segment_names, '[[', 8)
    
    segmentation_objects <- lapply(segmentation_files, function(x) get(load(x)))
    names(segmentation_objects) <- c(segment_names)  
  }  
}

load_bin_objs <- function(segmentation_files){
  if (length(segmentation_files) == 1){
    bin_objects <- get(load(segmentation_files))
    bod <- bin_objects$data
    bod <- bod[, !grepl("none|N_1.*N_1", colnames(bod))]
    bod <- bod %>% 
      dplyr::mutate(n_pos = rowSums(. > 0)) %>% 
      dplyr::mutate(n_neg = rowSums(. < 0)) %>%
      dplyr::select(chrom, maploc, n_pos, n_neg) %>% 
      identity()
    return(bod)
  } else{   # if loading rdata files containing single sample per file ---------------
    segment_names <- strsplit(segmentation_files,'/')
    segment_names <- sapply(segment_names, '[[', 8)
    
    bod <- lapply(segmentation_files, function(x) get(load(x)))
    names(bod) <- c(segment_names)  
  }  
}

load_all_objs <- function(segmentation_files){
  if (length(segmentation_files) == 1){
    SCNA_objects <- get(load(segmentation_files))
    return(SCNA_objects)
  } else{   # if loading rdata files containing single sample per file ---------------
    segment_names <- strsplit(segmentation_files,'/')
    segment_names <- sapply(segment_names, '[[', 8)
    
    SCNA_objects <- lapply(segmentation_files, function(x) get(load(x)))
    names(SCNA_objects) <- c(segment_names)  
  }  
}


RbindSegPlot <- function(SCNA_obj_list, sample_type = NULL, ...){
  SingleSegPlot <- function(seg_output, range.CNA = c(-2, 2),
                            color.palette = colorRampPalette(c("blue", "white", "red"))) {
    
    
    ## Use only tumor sample_ids
    if(sample_type == "tumor"){
      sample_ids <- unique(grep("^.*T.*N.*$", seg_output$ID, value = TRUE))
    } else if(sample_type == "cell_line"){   ## Use only cell line sample_ids
      sample_ids <- unique(grep("^.*CL.*N.*$", seg_output$ID, value = TRUE))
    } else
      # sample_ids <- c(unique(grep("^.*T.*N.*$", seg_output$ID, value = TRUE)), unique(grep("^.*CL.*N.*$", seg_output$ID, value = TRUE)))
      
      ## Use all sample_ids by default
      if(!exists("sample_ids")) {
        sample_ids <- unique(seg_output$ID)
      }
    
    ## Select sample_ids
    seg_output <- seg_output[seg_output$ID %in% sample_ids, ]
    
    # dna_copy_object$data <- dna_copy_object$data[, c("chrom", "maploc", sample_ids)]
    # names(dna_copy_object$data) <- c("chrom", "maploc", "seg.mean")
    
    order.sample_ids <- unique(seg_output$ID)
    seg_output$ID <- factor(seg_output$ID, levels = order.sample_ids)
    # return(list(seg_output, dna_copy_object$data))
    return(seg_output)
    
  }
  
  SCNA_obj_list <- lapply(SCNA_obj_list, SingleSegPlot) 
  # seg_obj_list <- lapply(SCNA_obj_list, "[[", 1)
  # point_obj_list <- lapply(SCNA_obj_list, "[[", 2)
  seg_obj_df <- rbindlist(SCNA_obj_list, idcol = "sample_id")
  # point_obj_df <- rbindlist(point_obj_list, idcol = "sample_id")
  # return(list("segments" = seg_obj_df, "bins" = point_obj_df))
  return(list("segments" = seg_obj_df))
}

create_seg_granges <- function(karyo_seg){
  
  chroms <- unique(karyo_seg$chrom)
  
  
  # prep data for karyoplots ----------------------------------------
  
  # rescale max gain of each sample to shrink small amp areas and increase dynamic range of gains 
  ceiling_gain <- 2
  total_max_gain <- max(karyo_seg$seg.mean)
  avg_max_gain <- group_by(karyo_seg, sample_id) %>% 
    summarize(top = max(seg.mean))
  avg_max_gain <- median(avg_max_gain$top)
  
  floor_loss <- -2
  total_min_loss <- min(karyo_seg$seg.mean)
  avg_min_loss <- group_by(karyo_seg, sample_id) %>% 
    summarize(top = min(seg.mean))
  avg_min_loss <- median(avg_min_loss$top)
  
  karyo_seg_gains <- karyo_seg %>% 
    filter(seg.mean > 0) %>% 
    mutate(seg.mean = ifelse(seg.mean > ceiling_gain, ceiling_gain, seg.mean)) %>%
    mutate(seg.mean = rescale(seg.mean, to=c(0,1)))
  
  karyo_seg_losses <- karyo_seg %>% 
    filter(seg.mean <= 0) %>% 
    mutate(seg.mean = ifelse(seg.mean < floor_loss, floor_loss, seg.mean)) %>%
    mutate(seg.mean = rescale(seg.mean, to=c(-1,0)))
  
  karyo_seg <- bind_rows(karyo_seg_gains, karyo_seg_losses)
  
  # karyo_seg <- mutate(karyo_seg, seg.mean = ifelse(seg.mean > avg_max_gain, avg_max_gain, seg.mean)) %>% 
  #   mutate(seg.mean = rescale(seg.mean, to=c(-1,1))) %>%  # rescale gains and losses to unit scale
  #   # mutate(seg.mean = ifelse(seg.mean < 0, rescale(seg.mean, to=c(-1,0)), seg.mean)) %>%  # rescale gains and losses to unit scale 
  #   # mutate(seg.mean = ifelse(seg.mean > 0, rescale(seg.mean, to=c(0,1)), seg.mean)) %>%  # rescale gains and losses to unit scale 
  #   identity()
  
  seg_obj_list <- split(karyo_seg, karyo_seg$sample_id) 
  seg_granges <- map(seg_obj_list, makeGRangesFromDataFrame, keep.extra.columns = TRUE)
  
  seg_names <- split(names(seg_granges), gsub("-.*", "", names(seg_granges)))
  seg_names <- lapply(seg_names, sort, decreasing = TRUE)
  seg_granges <- seg_granges[unlist(seg_names)]
  
  seg_granges <- seg_granges[!grepl("N", names(seg_granges))]
  
  
}


plot_base_seg <- function(cov_granges, tn, kp, chr_select, ...){
  
  # kp <- plotKaryotype(genome="hg19", plot.type = 1, plot.params = plot.params, labels.plotter = NULL, chromosomes = chr_select) %>%
  #   kpAddChromosomeNames(col="red", cex = 0.3) %>%
  #   kpAddBaseNumbers()
  
  
  tr.i <- 0.037
  tr.o <- 0.040
  
  print(unique(cov_granges$sample_id))
  kp <- kpDataBackground(kp, r0=tr.o*tn, r1=tr.o*tn+tr.i, data.panel = 2) %>%
    # kpAxis(ymin =0, y = 1, cex = 0.3, r0=(tr.o*tn), r1=(tr.o*tn+tr.i)) %>% 
    kpHeatmap(cov_granges, y = cov_granges$seg.mean, ymin=-1, ymax=1, r0=tr.o*tn, r1=tr.o*tn+tr.i, colors = c("blue", "white", "red"), data.panel = 2) %>% 
    kpAddLabels(labels=unique(cov_granges$sample_id), r0=tr.o*tn, r1=tr.o*tn+tr.i,  pos=1, label.margin = 0.04, col="red", cex=0.5, data.panel = 2) %>% 
    identity()
  
}

plot_all_to_file <- function(raw_cov_list, file_name, chr_select, marker_granges, ...) {
  tr.i <- 0.037
  tr.o <- 0.040
  num_seg_granges <- seq(0, length(raw_cov_list)-1)
  tn = max(num_seg_granges)+1
  png(file_name, width = 5, height = 4, units = 'in', res = 800)
  plot.params <- getDefaultPlotParams(plot.type=3)
  plot.params$ideogramheight <- 3
  plot.params$topmargin <- 5
  plot.params$data1height <- 1
  plot.params$data1inmargin <- 1
  plot.params$data2inmargin <- 1
  plot.params$data2height <- 150
  plot.params$data1outmargin <- 1
  plot.params$bottommargin <- 1
  plot.params$leftmargin <- 0.1
  kp <- plotKaryotype(genome="hg19", plot.type = 3, plot.params = plot.params, labels.plotter = NULL, chromosomes = chr_select) %>% 
    kpAddChromosomeNames(col="red", cex = 0.45, srt = 45) 
  map2(raw_cov_list, num_seg_granges, plot_base_seg, kp, chr_select, marker_granges)
  kpPlotRegions(kp, data = marker_granges, r0=tr.o*tn, r1=tr.o*tn+tr.i, col = "black", lty=1, lwd=0.25, border="black", data.panel=2) %>% 
    kpAddLabels(labels="kooi regions", r0=tr.o*tn, r1=tr.o*tn+tr.i,  pos=1, label.margin = 0.04, col="red", cex=0.5, data.panel = 2)
  dev.off()
}

load_seg_files <- function(seg_file){
  list.files(seg_file, "^segment.Rdata$", recursive = TRUE, full.names = TRUE)
}

seg_span_gene <- function(baf_segment, chr, feature_start, feature_end){
  baf_segment <- dplyr::filter(baf_segment, chrom == chr & (!(start < feature_start & end < feature_end) & !(start > feature_end & end > feature_end)))
  # baf_segment <- baf_segment[(baf_segment$chrom == chr & (!(baf_segment$start < feature_start & baf_segment$end < feature_end) & !(baf_segment$start > feature_end & baf_segment$end > feature_end)))]
}

clean_bin_objs <- function(segmentation_files, ...){
  bin_objects <- lapply(segmentation_files, load_bin_objs)
  bin_objects <- reduce(bin_objects, left_join, by = c("chrom", "maploc")) %>% 
    dplyr::mutate(n_pos = n_pos.x + n_pos.y) %>% 
    dplyr::mutate(n_neg = n_neg.x + n_neg.y) %>%
    dplyr::mutate(gain_loss_diff = n_pos - n_neg) %>%
    # dplyr::select(chrom, maploc, gain_loss_diff) %>%
    identity()
                                    
  return(bin_objects)
}

plot_bin_objects <- function(bin_objects, filename){
  ggplot(bin_objects, aes(x=maploc, y=gain_loss_diff, color=chrom)) +
    geom_line()
}

format_SCNA_data <- function(dna_copy_object, samples, remove_pattern, sample_type = NULL) {
  
  if(sample_type=="proband"){
    samples <- unique(grep("[[:digit:]]_1", dna_copy_object$output$ID, value = TRUE ))
  } else if (sample_type == "all"){
    samples <- unique(dna_copy_object$output$ID)
  }
  
  # none_refs <- unique(grep("none", samples, value = TRUE))
  # normal_normal <- unique(grep("N_1.*N_1", samples, value = TRUE))
  
  remove_samples <- unlist(lapply(remove_pattern, function(x) unique(grep(x, samples, value = TRUE))))
  ## Use all samples by default
  if(missing(samples)) {
    samples <- unique(dna_copy_object$output$ID)
  }
  
  ## Select samples
  dna_copy_object$output <- dna_copy_object$output[!dna_copy_object$output$ID %in% remove_samples, ]
  
  order.samples <- unique(dna_copy_object$output$ID)
  dna_copy_object$output$ID <- factor(dna_copy_object$output$ID, levels = order.samples)
  return(dna_copy_object$output)
  
}

find_seg_genes <- function(my_data, outfile = NULL) {
  my_data <- my_data %>% 
    dplyr::rename(sample = ID, start = loc.start, end = loc.end, seg.mean = seg.mean, chrom = chrom)
  
  peak <- my_data %>% 
    dplyr::filter(!is.na(seg.mean))
  
  peak <- peak %>% 
    as.data.frame %>%
    dplyr::select(chrom, start, end, sample, seg.mean) %>%
    mutate(chrom=paste0('chr', chrom)) %>%
    mutate(chrom=case_when(
      chrom == "chr23" ~ "chrX",
      chrom == "chr24" ~ "chrY",
      TRUE ~ as.character(chrom)
    )) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>% 
    identity()
  
  peak2 <- mergeByOverlaps(genes(txdb), peak)

  gene_overlaps <- data.frame(peak2$`genes(txdb)`)[,c("seqnames", "start", "end")]
  peaks <- data.frame(peak2$peak)[,c("start", "end")]
  colnames(peaks) <- c("seg.start", "seg.end")
  
  peaks <- cbind(gene_overlaps, peaks, "gene_id" = peak2$gene_id, "sample" = peak2$sample, "seg.mean" = peak2$seg.mean)
  
  entrez_to_symbol <- as.data.frame(org.Hs.egSYMBOL)
  
  peaks <- dplyr::right_join(entrez_to_symbol, peaks, by = "gene_id")
  
  # write.table(peaks, outfile, sep = ",", row.names = FALSE)
  return(peaks)
  
}

clean_gene_objs <- function(peak_genes_dat, ...){
  
  gene_objects <- group_by(peak_genes_dat, gene_id, seqnames, symbol, start, end) %>%  
    dplyr::count(sign(seg.mean)) %>% 
    tidyr::spread(`sign(seg.mean)`, n) %>% 
    dplyr::rename(n_neg = `-1`, n_pos = `1`) %>%
    dplyr::mutate_all(funs(replace(., is.na(.), 0))) %>%
    dplyr::mutate(gain_loss_diff = n_pos - n_neg) %>%
    dplyr::select(gene_id, symbol, seqnames, start, end, gain_loss_diff) %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::filter(seqnames %in% c("chr1", "chr2", "chr6", "chr7", "chr13", "chr16", "chr19")) %>%
    identity()
  
  return(gene_objects)
}

plot_gene_objects <- function(gene_objects, title, filename){

  gene_objects <- clean_gene_objs(gene_objects)
  gene_plot <- ggplot(gene_objects, aes(x=start, y=gain_loss_diff, color=seqnames)) +
    geom_line() +
    ggtitle(title)
  
  print(gene_plot)
  return(gene_plot)
  # 
  # ggplot(gene_objects, aes(x=start, y=gain_loss_diff, color=seqnames)) +
  #   geom_segment(aes(x=start, xend=end, y=gain_loss_diff, yend=gain_loss_diff)) +
  #   ggtitle(title)
  
}


# load segmentation data ---------------------------------------------------------------

# load kooi SCNA peak regions 
kooi_peak_regions <- read.table("~/rb_pipeline/doc/SCNA/kooi_SCNA_peak_regions.csv", stringsAsFactors = F)

kooi_peak_granges <- makeGRangesFromDataFrame(kooi_peak_regions)

if (cohort == "vc"){
  old_seg_files = "~/rb_pipeline/output/copywriter/100kb/old_CNAprofiles/"
  new_seg_files = "~/rb_pipeline/output/copywriter/100kb/new_CNAprofiles/"
  segmentation_files <- c(old_seg_files, new_seg_files)
  segmentation_files <- lapply(segmentation_files, load_seg_files)
  segmentation_objects <- do.call(c, lapply(segmentation_files, load_seg_objs))
  
  # tidy segmentation data ---------------------------------------------------------------
  multi_seg_out_tumor_stchlk <- RbindSegPlot(segmentation_objects, sample_type = "tumor")
  # multi_seg_out_cell_line_stchlk <- RbindSegPlot(segmentation_cl, sample_type = "cell_line")
  multi_seg_out_all_stchlk <- RbindSegPlot(segmentation_objects, sample_type = "all")
  
  karyo_seg <- multi_seg_out_all_stchlk$segments %>% 
    dplyr::rename("start" = "loc.start", "end" = "loc.end") %>% 
    mutate(chrom = paste0("chr", chrom)) %>% 
    mutate(chrom = gsub("chr23", "chrX", chrom)) %>% 
    mutate(chrom = gsub("chr24", "chrY", chrom)) %>% 
    mutate(sample_id = gsub("log2.", "", gsub("_.*", "", ID))) %>% 
    identity()
  
  # karyo_bin <- multi_seg_out_all_stchlk$bins %>% 
  #   dplyr::rename("start" = "loc.start", "end" = "loc.end") %>% 
  #   mutate(chrom = paste0("chr", chrom)) %>% 
  #   mutate(chrom = gsub("chr23", "chrX", chrom)) %>% 
  #   mutate(chrom = gsub("chr24", "chrY", chrom)) %>% 
  #   mutate(sample_id = gsub("log2.", "", gsub("_.*", "", ID))) %>% 
  #   identity()
  
  karyo_seg <- karyo_seg[!grepl("none", karyo_seg$ID),]
  karyo_seg <- karyo_seg[!grepl("N", karyo_seg$sample_id),]
  
  seg_granges <- create_seg_granges(karyo_seg)
  
  seg_granges <- sort_T_b4_CL(seg_granges)
  

  # plot vc seg
  plot_all_to_file(seg_granges, "~/rb_pipeline/doc/SCNA/chrom_vc_seg_coverage.png", chr_select = "auto", marker_granges = kooi_peak_granges)
  
  # output vc seg
  vc_seg_path <- "~/rb_pipeline/results/SCNA/vc_SCNA_segments.txt"
  
  write.table(karyo_seg, vc_seg_path, sep = "\t")
  
} else if (cohort == "reynolds"){
  
  reynolds_seg_files = "~/rb_pipeline/output/copywriter/20kb/CNAprofiles/" 
  segmentation_files <- reynolds_seg_files
  segmentation_files <- lapply(segmentation_files, load_seg_files)
  segmentation_objects <- do.call(c, lapply(segmentation_files, load_seg_objs))
  
  # tidy segmentation data ---------------------------------------------------------------
  multi_seg_out_all_stchlk <- RbindSegPlot(segmentation_objects, sample_type = "all")
  
  karyo_seg <- multi_seg_out_all_stchlk$segments %>% 
    dplyr::rename("start" = "loc.start", "end" = "loc.end") %>% 
    mutate(chrom = paste0("chr", chrom)) %>% 
    mutate(chrom = gsub("chr23", "chrX", chrom)) %>% 
    mutate(chrom = gsub("chr24", "chrY", chrom)) %>% 
    mutate(sample_id = gsub("log2.", "", gsub("_.*", "", ID))) %>% 
    identity()
  
  karyo_seg <- karyo_seg[grepl("none", karyo_seg$ID),]
    
  seg_granges <- create_seg_granges(karyo_seg)

  # plot reynolds seg
  plot_all_to_file(seg_granges, "~/rb_pipeline/doc/SCNA/chrom_reynolds_seg_coverage.png", chr_select = "auto", marker_granges = kooi_peak_granges)

  # output reynolds seg
  reynolds_seg_path <- "~/rb_pipeline/results/SCNA/reynolds_SCNA_segments.txt"
  write.table(karyo_seg, reynolds_seg_path, sep = "\t")
  
}


# plot gain loss difference in exome samples ------------------------------


if (cohort == "vc"){
  # load SCNA seg data
  SCNA_objects <- lapply(segmentation_files, load_all_objs)
  SCNA_data <- purrr::map_df(SCNA_objects, format_SCNA_data, remove_pattern = c("none", "N_1.*N_1"), sample_type="all")
  chroms <- unique(SCNA_data$chrom)
  
  # find genes present in SCNA segments
  peak_genes_dat <- find_seg_genes(SCNA_data)
  
  # filter peak_genes by sample type tumor or cell line
  T_peak_dat <- peak_genes_dat[grepl("T", peak_genes_dat$sample),]
  CL_peak_dat <- peak_genes_dat[grepl("CL", peak_genes_dat$sample),]
  
  # plot gain loss differences in sample types
  plot_gene_objects(T_peak_dat, "SCNA peak regions in VC Tumors")
  ggsave("~/rb_pipeline/doc/SCNA/SCNA_meta_plot_vc_tumor.png")
  plot_gene_objects(CL_peak_dat, "SCNA peak regions in VC Cell Lines")
  ggsave("~/rb_pipeline/doc/SCNA/SCNA_meta_plot_vc_cell_line.png")
  
} else if (cohort == "reynolds"){
  # load SCNA seg data
  SCNA_objects <- lapply(segmentation_files, load_all_objs)
  SCNA_data <- purrr::map_df(SCNA_objects, format_SCNA_data, remove_pattern = c("CL_1.*CL_1"), sample_type="all")
  chroms <- unique(SCNA_data$chrom)
  
  # find genes present in SCNA segments
  peak_genes_dat <- find_seg_genes(SCNA_data)
  
  # filter peak_genes by sample type tumor or cell line
  CL_peak_dat <- peak_genes_dat[grepl("CL", peak_genes_dat$sample),]
  
  # plot gain loss differences in sample types
  plot_gene_objects(CL_peak_dat, "SCNA peak regions in Reynolds Cell Lines")
  ggsave("~/rb_pipeline/doc/SCNA/SCNA_meta_plot_reynolds_cell_line.png")
}

# find segments spanning genes of interest --------------------------------

rb1_start = 48877883
rb1_end = 49056026

mycn_start <- 16080683
mycn_end <- 16087129

rb_segments <- seg_span_gene(karyo_seg, "chr13", rb1_start, rb1_end) %>% 
  filter(seg.mean < -0.5)

mycn_segments <- baf_rb <- seg_span_gene(karyo_seg, "chr2", mycn_start, mycn_end) %>% 
  filter(seg.mean > 1)

# segmentation_files = list.files("~/rb_pipeline/output/copywriter/100kb/CNAprofiles/", "^segment.Rdata$", recursive = TRUE, full.names = TRUE)

# if loading an rdata file containing multiple samples -------------------------

# saveRDS(multi_seg_out_all_stchlk, "~/rb_pipeline/results/SCNA/multi_seg_out_reynolds.rds")

# karyo_bins <- multi_seg_out_all_stchlk$bins %>% 
#   dplyr::rename("start" = "maploc") %>%
#   mutate(start  = as.numeric(start)) %>% 
#   mutate(end = (start + 20000)) %>% 
#   mutate(chrom = paste0("chr", chrom)) %>% 
#   mutate(chrom = gsub("chr23", "chrX", chrom)) %>% 
#   mutate(chrom = gsub("chr24", "chrY", chrom)) %>%
#   mutate(seg.mean = as.numeric(seg.mean)) %>% 
#   mutate(seg.mean = rescale(seg.mean, to = c(-1,1))) %>% 
#   filter(!is.na(seg.mean))

# bin_obj_list <- split(karyo_bins, as.factor(karyo_bins$sample_id)) 
# bin_granges <- map(bin_obj_list, makeGRangesFromDataFrame, keep.extra.columns = TRUE)
# 
# bin_names <- split(names(bin_granges), gsub("-.*", "", names(bin_granges)))
# bin_names <- lapply(bin_names, sort, decreasing = TRUE)
# bin_granges <- bin_granges[unlist(bin_names)]
# saveRDS(bin_granges, "~/rb_pipeline/results/SCNA/bin_granges.rds")


# prep functions for karyoplot --------------------------------------------

plot_path = "~/rb_pipeline/results/SCNA/karyoploter/chrom_16_coverage.png"

# plot_base_seg(seg_granges$`14-CL`, 1, kp, "auto", markers)


# plot all-chromosome karyoplots ------------------------------------------

# plot_all_to_file(bin_granges, "~/rb_pipeline/doc/SCNA/chrom_all_bin_coverage.png", chr_select = "auto", markers)


# plot_all_to_file(bin_granges, "./chrom_all_bin_coverage.png", chr_select = "auto", markers)
# plot_all_to_file(seg_granges, "./chrom_all_seg_coverage.png", chr_select = "auto", markers)



# plot chromosome-specific karyoplots (default chr 16) --------------------
markers <- data.frame(chr="chr16", pos=(53468351), labels="RBL2")
RBL2.region <- toGRanges(data.frame("chr16", 53050000, 53650000))

plot_chrom_to_file <- function(raw_cov_list, file_name, chr_select, marker_df, zoom_region = NULL, tick_dist = 10000000, ...) {
  
  num_seg_granges <- seq(0, length(raw_cov_list)-1)
  png(file_name, width = 4, height = 4, units = 'in', res = 800)
  plot.params <- getDefaultPlotParams(plot.type=3)
  plot.params$ideogramheight <- 5
  plot.params$data1height <- 20
  plot.params$data1inmargin <- 1
  plot.params$data2height <- 150
  plot.params$ttommargin <- 20
  plot.params$topmargin <- 20
  kp <- plotKaryotype(genome="hg19", plot.type = 3, plot.params = plot.params, labels.plotter = NULL, chromosomes = chr_select, zoom = zoom_region) %>% 
    kpAddChromosomeNames(col="red", cex = 0.3) %>% 
    kpPlotMarkers(chr=markers$chr, x=markers$pos, y = 0.6, labels=markers$labels, data.panel = 1, cex = 0.3) %>% 
    kpAddBaseNumbers(tick.dist = tick_dist, minor.tick.dist = tick_dist)
  
  map2(raw_cov_list, num_seg_granges, plot_base_seg, kp, chr_select, marker_df)
  dev.off()
}

plot_chrom_to_file(seg_granges, "~/rb_pipeline/doc/SCNA/karyoploter/chrom_16_seg_coverage.png", chr_select = "chr16", markers, zoom_region = NULL)

chrom_16_bin_granges <- lapply(bin_granges, function(x){x[seqnames(x) == "chr16"]})
RBL2_granges <- lapply(bin_granges, subsetByOverlaps, RBL2.region, type = "within")

plot_chrom_to_file(chrom_16_bin_granges, "~/rb_pipeline/doc/SCNA/karyoploter/chrom_16_bin_coverage.png", chr_select = "chr16", markers, zoom_region = NULL)
plot_chrom_to_file(RBL2_granges, "~/rb_pipeline/doc/SCNA/karyoploter/RBL2_bin_coverage.png", chr_select = "chr16", markers, zoom_region = RBL2.region, tick_dist = 80000)

# biomaRt look up genes ---------------------------------------------------

chrom_16_loss_regions <- lapply(bin_granges, function(x){x[seqnames(x) == "chr16"]})
chrom_16_loss_regions <- lapply(chrom_16_bin_granges, function(x){seqlevelsStyle(x) <- "Ensembl"; return(x)})

chrom_all_loss_regions <- lapply(bin_granges, function(x){seqlevelsStyle(x) <- "Ensembl"; return(x)})

lookup_genes <- function(my_data, threshold = NULL) {
  
  gene_sym_from_coord <- function(my_data){
    
    my_data <- my_data[(elementMetadata(my_data)[, "seg.mean"] < threshold[1] | elementMetadata(my_data)[, "seg.mean"] > threshold[2])]
    if(length(my_data) > 0){
      results <- getBM(
        attributes = c("hgnc_symbol", "chromosome_name", "start_position","end_position"),
        filters = c("chromosomal_region"),
        values=my_data, 
        mart = ensembl) 
      
    }
  }
  
  peak_genes <- invisible(lapply(my_data, gene_sym_from_coord))
  peak_genes[sapply(peak_genes, is.null)] <- NULL
  # peak_genes <- rbindlist(peak_genes)
  # 
  # peak_genes <- peak_genes %>% 
  #   filter(!(start_position < start & end_position < start) & !(start_position > end & end_position > end)) %>% 
  #   arrange(seg.mean)
  
  return(peak_genes)
  
}

chr16_genes <- lookup_genes(chrom_16_loss_regions, c(-0.5, 0.5))
chr_all_genes <- lookup_genes(chrom_all_loss_regions, c(-0.99,0.99))

tidy_lookup_genes <- function(gene_set, csv_out){
  genes_df <- rbindlist(gene_set, idcol = "sample") %>% 
    filter(hgnc_symbol != "") %>% 
    group_by(hgnc_symbol) %>% 
    mutate(recurrence = paste(as.character(sample), collapse = ";")) %>% 
    mutate(counts = n()) %>% 
    filter(row_number() == 1)
  write.table(genes_df, csv_out, sep = "\t", row.names = FALSE) 
  return(genes_df)
}

chr16_peak_genes <- tidy_lookup_genes(chr16_genes, "~/rb_pipeline/results/SCNA/chr16_peak_genes.csv")

chr_all_peak_genes <- tidy_lookup_genes(chr_all_genes, "~/rb_pipeline/results/SCNA/chr_all_peak_genes.csv")




