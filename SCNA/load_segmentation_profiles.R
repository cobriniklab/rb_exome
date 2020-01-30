#!/usr/bin/Rscript

library(reshape2)


segmentation_files = list.files("hi_res_copywriter", "^segment.Rdata$", recursive = TRUE, full.names = TRUE)

segment_names       <- strsplit(segmentation_files,'/')
segment_names   <- unlist(segment_names)[4*(1:length(segmentation_files))-2]

segmentation_objects <- lapply(segmentation_files, function(x) get(load(x)))
names(segmentation_objects) <- c(segment_names)

OverviewPlot <- function(DNAcopy.object, samples, range.CNA = c(-2, 2),
                         color.palette = colorRampPalette(c("blue", "white", "red"))(49)) {
  
    
  p1 = paste(as.character(DNAcopy.object), "_plot")
  
    ## Use only tumor samples
    if(samples == "tumor"){
      samples <- unique(grep("^.*T.*N.*$", DNAcopy.object$output$ID, value = TRUE))
      print(samples)
    }
    
    ## Use only cell line samples
    if(samples == "cell_line"){
      samples <- unique(grep("^.*CL.*N.*$", DNAcopy.object$output$ID, value = TRUE))
      print(samples)
    }
    
    ## Select samples
    DNAcopy.object$output <- DNAcopy.object$output[DNAcopy.object$output$ID %in% samples, ]
    
    ## Cap range
    DNAcopy.object$output$seg.mean[DNAcopy.object$output$seg.mean < range.CNA[1]] <- range.CNA[1]
    DNAcopy.object$output$seg.mean[DNAcopy.object$output$seg.mean > range.CNA[2]] <- range.CNA[2]
    
    ## Reshape data.frame according to sample name
    # reshape2::dcast -> to wider data.frame; right-hand of tilde needs to be a 'measured variable' (i.e., needs to go into columns)
    # reshape2::melt -> to narrower data.frame
    # Names are changed by dcast according to level order -> change order levels (!)
    order.samples <- unique(DNAcopy.object$output$ID)
    DNAcopy.object$output$ID <- factor(DNAcopy.object$output$ID, levels = order.samples)
    DNAcopy.object.cast <- reshape2::dcast(data = DNAcopy.object$output, formula = ID + chrom + loc.start + loc.end + num.mark + seg.mean ~ ID
                                           , value.var = "num.mark")
    
    DNAcopy.object.cast[is.na(DNAcopy.object.cast)] <- 0
    
    ## Calculate number of samples
    sample.number <- ncol(DNAcopy.object.cast) - 6
    
    ## Collapse data
    DNAcopy.object.cast.aggregate <- aggregate(DNAcopy.object.cast[, c("chrom", "num.mark")], by = list(DNAcopy.object.cast$chrom), FUN = sum)
    
    ## Calculate scaling factors
    range.factor <- range.CNA[2] - range.CNA[1]
    
    ## Create overviewPlot
    par(mfrow = c(1, 2 + sample.number), mar = c(7, 0, 1, 0))
    barplot(matrix(DNAcopy.object.cast.aggregate[, "num.mark"], dimnames = list(NULL, "Chr")), col = c("black", "white"), border = NA
            , yaxt = "n", xlim = c(0, 0.1), width = 0.1, las = 2)
    
    for(i in 1:sample.number) {
      # color_scale = colorRampPalette(1 + (DNAcopy.object.cast$seg.mean - range.CNA[1]) / range.factor * (length(color.palette) - 1))
      #  p1 <- ggplot(DNAcopy.object.cast[, 6 + i, drop = FALSE], aes(sample)) + geom_bar(fill = color_scale)    
     # print(p1)
      barplot(as.matrix(DNAcopy.object.cast[, 6 + i, drop = FALSE])
             , col = color.palette[1 + (DNAcopy.object.cast$seg.mean - range.CNA[1]) / range.factor * (length(color.palette) - 1)], border = NA
              , yaxt = "n", xlim = c(0, 0.1), width = 0.1, las = 2)
    }
   
    #barplot(matrix(rep(1, length(color.palette)), ncol = 1, dimnames = list(NULL, paste(range.CNA, collapse = " to ")))
    #      , col = color.palette, border = NA, yaxt = "n", beside = FALSE, xlim = c(0, 0.1), width = 0.1, las = 2)
  
}

OverviewPlot(segmentation_objects$`29-CL_1`, "cell_line", range.CNA = c(-2, 2))
OverviewPlot(segmentation_objects$`49-T_1`, "tumor", range.CNA = c(-2, 2))
OverviewPlot(segmentation_objects, "cell_line", range.CNA = c(-2, 2))
plots = list()
plots <- lapply(segmentation_objects, function(x){ OverviewPlot(x, "cell_line", range.CNA = c(-2,2))})









