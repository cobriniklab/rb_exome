#' Make a waterfall plot
#'
#' @param my_variants 
#' @param clindata 
#' @param gene_lab_size 
#'
#' @return
#' @export
#'
#' @examples
make_waterfall_plot <- function(my_variants, clindata, subset_genes = NULL, sample_order = NULL, gene_lab_size = 10, mainLabelSize = 4){
  # in long format with colnames = c("sample", "variable", "value")
  
  if (!is.null(subset_genes)){
    my_variants <- dplyr::filter(my_variants, gene %in% subset_genes)
  }
  
  clindata <- clindata %>% 
    tidyr::gather("variable", "value", -sample) %>% 
    dplyr::filter(!grepl("age|author", variable)) %>% 
    dplyr::filter(sample %in% unique(my_variants$sample)) %>% 
    # rename(sample = Sample) %>%
    identity()
  
  clinvars <- unique(clindata$value)
  
  clinVarCol <- paletteer::paletteer_d("ggsci", "default_igv")[1:length(clinvars)] %>% 
    set_names(clinvars)
  
  clinVarOrder <- clinvars 
  
  input_wat <- my_variants %>% 
    ungroup() %>% 
    dplyr::select(sample, gene, 
                  variant_class = Consequence, 
                  amino.acid.change = HGVSp) %>% 
    dplyr::mutate(variant_class = replace_na(variant_class, "unknown")) %>% 
    dplyr::mutate(sample = factor(sample, levels = sort(unique(.$sample)))) %>% 
    identity()
  
  mutation_priority <- as.character(unique(input_wat$variant_class))
  grob1 <- waterfall(input_wat,
                     fileType = "Custom",
                     variant_class_order = mutation_priority, 
                     out = "grob",
                     mainPalette = paletteer::paletteer_d("ggsci", "default_igv")[1:length(unique(input_wat$variant_class))],
                     clinData = clindata,
                     clinLegCol = 5,
                     clinVarCol = clinVarCol,
                     clinVarOrder = clinVarOrder,
                     mainXlabel = TRUE,
                     main_geneLabSize = gene_lab_size,
                     mainLabelSize = mainLabelSize,
                     plotMutBurden = FALSE,
                     plot_proportions = FALSE, 
                     sampOrder = sample_order
  )
  
  return(grob1)
  
}

#' prep stachelek clindata
#'
#' @param clindata 
#'
#' @return
#' @export
#'
#' @examples
prep_stachelek_clindata <- function(clindata){
  clindata[is.na(clindata)] <- 0
  
  # clindata_cl = mutate(clindata, sample = paste0(sample, "-CL")) %>% 
  #   identity()
  # 
  # clindata_t = mutate(clindata, sample = paste0(sample, "-T")) %>% 
  #   identity()
  # 
  # clindata = dplyr::bind_rows(clindata_cl, clindata_t)
  
  clindata <- 
    clindata %>% 
    dplyr::mutate(sample_suffix = case_when(series == "CHLA-VC-RB" ~ "CL;T",
                                          series == "CHLA-RB" ~ "CL")) %>%
    dplyr::mutate(sample = gsub(".*RB-|CHLA-", "", sample)) %>%
    dplyr::mutate(sample_suffix = stringr::str_split(sample_suffix, ";")) %>%
    tidyr::unnest(sample_suffix) %>%
    dplyr::mutate(sample = paste0(sample, "-", sample_suffix)) %>% 
    # dplyr::mutate(sample = gsub(".*RB-|CHLA-", "", sample)) %>% 
    identity()
  
}

extract_clindata <- function(variants, clindata) {
  clindata_cols <- c(colnames(clindata), "author")
  
  clindata = dplyr::ungroup(variants) %>% 
    dplyr::select(one_of(clindata_cols)) %>% 
    dplyr::mutate(series = dplyr::coalesce(series, author)) %>% 
    dplyr::distinct()
}

list_recurrent_genes <- function(variant_df){
  variant_df <- 
    prior_author_variants %>% 
    group_by(sample, gene) %>%
    dplyr::filter(row_number() == 1) %>% 
    group_by(gene) %>% 
    summarize(n = dplyr::n()) %>% 
    mutate(freq = n / number_of_samples) %>%
    # filter(gene == "BCOR") %>% 
    mutate(total = number_of_samples) %>% 
    dplyr::arrange(desc(n)) %>% 
    identity()
  
} 

#' Run WebGestaltR
#'
#' @param var_df 
#' @param gene_list_path 
#'
#' @return
#' @export
#'
#' @examples
run_webgestaltr <- function(var_df, gene_list_path) {
  gene_list_out <- fs::path(path_dir(gene_list_path), paste0(path_file(path_ext_remove(gene_list_path)), "_results"))
  dir_create(gene_list_out)
  
  write_delim(var_df["gene"], gene_list_path)
  
  test0 <- WebGestaltR(interestGeneFile = gene_list_path,
                       interestGeneType = "genesymbol",
                       referenceSet = "genome",
                       enrichDatabase = "geneontology_Biological_Process",
                       outputDirectory = gene_list_out,
                       sigMethod = "top",
                       topThr = 25)
}

filter_variant_set  <- function (my_variants, datatype, gnomad_threshold = 0.0005) {
  if (datatype %in% c("tn")) {
    # browser()
    myv <- filter_calls(my_variants)
    myv <- myv %>%
      dplyr::filter((AD.NORMAL.1 != 0 | AD.NORMAL.2 != 0)) %>% # filter out variants with zero normal reads (ref or alt)
      dplyr::filter(gnomad.AF < gnomad_threshold | is.na(gnomad.AF)) %>%  # filter out variants with a gnomad frequency greater than set threshold
      identity()
    
    return(myv)
    
  } else if (datatype == "pon") {
    # browser()
    myv <- filter_calls(my_variants)
    myv <- myv %>%
      dplyr::filter(gnomad.AF < gnomad_threshold | is.na(gnomad.AF)) %>%  # filter out variants with a gnomad frequency greater than set threshold
      identity()
    
    return(myv)
  } else if (datatype == "strelka") {
    # browser()
    myv <- filter_calls(my_variants)
    myv <- myv %>%
      dplyr::filter(gnomad.AF < gnomad_threshold | is.na(gnomad.AF)) %>%  # filter out variants with a gnomad frequency greater than set threshold
      identity()
    
    return(myv)
  } else {
    stop("'datatype' describes the method of variant calling, either tumor/normal paired (tn) or panel of normals (pon)")
  }
}
