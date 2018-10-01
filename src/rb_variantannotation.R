#!/usr/bin/Rsript


# load required libraries -------------------------------------------------


library(VariantAnnotation)
library(tidyverse)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(PolyPhen.Hsapiens.dbSNP131)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(COSMIC.67)
library(biobroom)
library(MafDb.gnomADex.r2.0.1.hs37d5)
library(biomaRt)
ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")

#clinvar <-  read.table('../Homo_sapiens/ClinVar/clinvar_alleles.single.b37.tsv',sep='\t',header=T,quote='',comment.char='')


# load input --------------------------------------------------------------


mutect2_tumor_filenames <- list.files(path="./output/gatk", pattern="[[:digit:]]-T_1_filtered.vcf.gz$", full.names = TRUE)
mutect2_cell_line_filenames <- list.files(path="./output/gatk", pattern="[[:digit:]]-CL_1_filtered.vcf.gz$", full.names = TRUE)
mutect2_pon_filenames <- list.files(path="./output/gatk", pattern="[[:digit:]]-N_1_pon_annotated.vcf.gz$", full.names = TRUE)
mutect2_reynolds_filenames <- list.files(path="./output/gatk", pattern="[[:digit:]]{3}-.*rb_anno.vcf.gz$", full.names = TRUE)


strelka_snv_filenames <- list.files(path="./output/strelka/", pattern="somatic.snvs.vcf.gz$",  full.names = TRUE, recursive = TRUE)
strelka_indel_filenames <- list.files(path="./output/strelka/", pattern="somatic.indels.vcf.gz$",  full.names = TRUE, recursive = TRUE)

mutect2_tumor_names <-substr(mutect2_tumor_filenames,15,18) 
mutect2_cell_line_names <-substr(mutect2_cell_line_filenames,15,19) 
mutect2_normal_names <-substr(mutect2_pon_filenames,15,18)
mutect2_reynolds_names <- substr(mutect2_reynolds_filenames,15,17)

strelka_snv_names <- substr(strelka_snv_filenames ,19,22)
strelka_indel_names <- substr(strelka_indel_filenames ,19,22)


# read vcf ----------------------------------------------------------------

#svp <- ScanVcfParam(info="AF", geno="GT")

mutect2_tumor_list <- lapply(mutect2_tumor_filenames, function(x)readVcf(x, "hg19"))
mutect2_cell_line_list <- lapply(mutect2_cell_line_filenames, function(x)readVcf(x, "hg19"))
mutect2_pon_list <- lapply(mutect2_pon_filenames, function(x)readVcf(x, "hg19"))
mutect2_reynolds_list <- lapply(mutect2_reynolds_filenames, function(x)readVcf(x, "hg19")) 

strelka_snv_list <- lapply(strelka_snv_filenames, function(x)readVcf(x, "hg19"))
strelka_indel_list <- lapply(strelka_indel_filenames, function(x)readVcf(x, "hg19"))

fl <- mutect2_tumor_filenames[[1]]

names(mutect2_tumor_list) <- mutect2_tumor_names
names(mutect2_cell_line_list) <- mutect2_cell_line_names
names(mutect2_pon_list) <- mutect2_normal_names
names(mutect2_reynolds_list) <- mutect2_reynolds_names


names(strelka_snv_list) <- strelka_snv_names
names(strelka_indel_list) <- strelka_indel_names

vcf <-  mutect2_tumor_list[[1]]

mutect2_tumor_pon_list <- c(mutect2_pon_list, mutect2_tumor_list)
mutect2_cell_line_pon_list <- c(mutect2_pon_list, mutect2_cell_line_list)

test_tumor_list <- lapply(mutect2_tumor_list, '[', c(1:10))
test_reynolds_list <- lapply(mutect2_reynolds_list, '[', c(1:10))

single_tidy <- function(my_vcf){

  evcf <- S4Vectors::expand(my_vcf)  
  df_all <- data.frame(rowRanges(evcf)) %>% 
      rownames_to_column("snp_id")
    

# predict coding/promoters ----------------------------------------------------------

  library(BSgenome.Hsapiens.UCSC.hg19)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  all_vars <- rowRanges(my_vcf)
  coding <- predictCoding(my_vcf, txdb, seqSource=Hsapiens)
  genes_txdb <- genes(txdb)
  promoters_txdb <- promoters(genes_txdb)
  promoters_anno <- subsetByOverlaps(coding, promoters_txdb)
  promoter <- tidy.GRanges(promoters_anno) %>% 
    dplyr::rename(seqnames = seqname) %>% 
    dplyr::mutate(region = "promoter") %>% 
    dplyr::select(seqnames, start, end, region)
  

# predict coding/frameshift -----------------------------------------------


  df_coding <- biobroom::tidy.GRanges(coding) %>%
    dplyr::select(seqname, start, end, GENEID, CONSEQUENCE, REFAA, VARAA) %>%
    dplyr::rename("seqnames" = "seqname") %>%
    left_join(df_all, by = c("seqnames", "start", "end")) %>%
    left_join(promoter, by = c("seqnames", "start", "end")) %>%
    mutate(region = ifelse(is.na(region), "coding", "promoter"))
  

# gnomad ------------------------------------------------------------------
  ncbiStyle <- mapSeqlevels(seqlevels(all_vars),"NCBI")
  all_vars <- renameSeqlevels(all_vars, ncbiStyle)
  all_vars <- all_vars[width(all_vars) == 1]
  mafdb <- MafDb.gnomADex.r2.0.1.hs37d5
  all_vars_vec <- unique(all_vars)
  mfdf <- mafByOverlaps(mafdb, all_vars_vec) 
  
  mfdf <-  biobroom::tidy.GRanges(mfdf)  
  mfdf <- mfdf %>% 
    dplyr::rename(seqnames = seqname) %>% 
    mutate(gnomad.AF = AF) %>% 
    mutate(seqnames = paste0("chr", seqnames))

  
  variants <- full_join(df_coding, mfdf, by = c("seqnames", "start", "end")) %>% 
    group_by(seqnames, start, end) %>% 
    filter(row_number() == 1) %>% 
    ungroup()
    
}


collate_vcfs <- function(vcf_list, flavor){
  browser()
  evcf_list <- lapply(vcf_list, function(x)single_tidy(x))
  if (flavor == "tn"){
    evcf_list0 <- evcf_list[1:10]
    evcf_list0 <- lapply(evcf_list0, function(x){
      colnames(x)[16:20] = c("GT.NORMAL", "AF.NORMAL", "AD.NORMAL.1", "AD.NORMAL.2", "DP.NORMAL")
      return(x)
    })
    evcf_list <- c(evcf_list, evcf_list[11:length(evcf_list)])
  } else if (flavor == "pon"){
    evcf_list <- evcf_list 
  } else {
    stop("'flavor' describes the method of variant calling, either tumor/normal paired (tn) or panel of normals (pon)")
  }
  
  
  tidy_vcfs <- data.table::rbindlist(evcf_list0, idcol = "sample", fill = TRUE) 
  
  # rfpred ------------------------------------------------------------------
  
  
  library(rfPred)
  
  rfp_input <- dplyr::select(tidy_vcfs, chrom = seqnames, pos = start, ref = REF, alt = ALT) %>% 
    mutate(chrom = gsub("chr", "", chrom))
  
  rfp0 <- rfPred_scores(variant_list=rfp_input,
                        data="~/rb_pipeline/bin/all_chr_rfPred.txtz",
                        index="~/rb_pipeline/bin/all_chr_rfPred.txtz.tbi")
  
  rfp0 <- mutate(rfp0, chromosome = paste0("chr", chromosome))
  
  tidy_vcfs <- dplyr::left_join(tidy_vcfs, rfp0, by = c("seqnames" = "chromosome", "start" = "position_hg19", "REF" = "reference", "ALT" = "alteration"))
  
}

tidy_tumors0 <- collate_vcfs(mutect2_tumor_pon_list, "tn")
tidy_cell_lines0 <- collate_vcfs(mutect2_cell_line_list, "tn")
tidy_normals0 <- collate_vcfs(mutect2_pon_list, "pon")
tidy_total0 <- collate_vcfs(mutect2_total_list, "tn")
tidy_reynolds0 <- collate_vcfs(mutect2_reynolds_list, "pon")


saveRDS(tidy_tumors0, "./results/tidy_mutect2_tumor_variants.rds")
saveRDS(tidy_normals0, "./results/tidy_mutect2_pon_variants.rds")
saveRDS(tidy_cell_lines0, "./results/tidy_mutect2_cell_line_variants.rds")
saveRDS(tidy_total0, "./results/tidy_mutect2_both_tumor_and_cell_line_variants.rds")
saveRDS(tidy_reynolds0, "./results/tidy_mutect2_reynolds_cell_line_variants.rds")

retidy_pon <- function(my_vcfs){
browser()
  my_vcfs <- dplyr::arrange(my_vcfs, desc(AF.NORMAL)) %>% 
    dplyr::group_by(sample, snp_id) %>% 
    dplyr::filter(row_number() == 1) %>%
    dplyr::group_by(snp_id) %>% 
    dplyr::mutate(recurrence=paste(sample,collapse=';')) %>%
    dplyr::mutate(counts = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(GENEID = as.character(GENEID)) 
  
  org_translation <- select(org.Hs.eg.db, as.character(my_vcfs$GENEID), c("SYMBOL", "GENENAME"))
  
  my_vcfs <- left_join(my_vcfs, org_translation, by = c("GENEID" = "ENTREZID"))
  
  variants <- my_vcfs %>% 
    group_by(snp_id) %>% 
    filter(row_number() == 1) %>% 
    ungroup() %>% 
    arrange(desc(counts))
  
  genes <- variants %>%
    ungroup() %>% 
    group_by(GENEID) %>%
    mutate(counts = sum(counts)) %>%
    mutate(recurrence = sapply(recurrence, c)) %>%
    mutate(recurrence = sapply(strsplit(as.character(recurrence) ,";"), function(x) paste(unique(x), collapse="; "))) %>% 
    mutate(gene_counts = sapply(strsplit(as.character(recurrence) ,";"), function(x) length(unique(x)))) %>%
    filter(row_number() == 1) %>%
    ungroup() %>% 
    arrange(desc(counts))
  
  
  samples <- variants %>% 
    ungroup() %>% 
    dplyr::select(sample, SYMBOL, GENENAME) %>% 
    dplyr::arrange(desc(sample)) 
  
  
  samples <- aggregate(SYMBOL ~ sample, data = unique(samples), paste, collapse = ",")
  
  
  
  
  
  my_list <-  list("variants" = variants, "genes" = genes, "samples" = samples)
  
}

retidy_vcfs <- function(my_vcfs){

  my_vcfs <- arrange(my_vcfs, desc(AF.TUMOR)) %>% 
    group_by(sample, snp_id) %>% 
    filter(row_number() == 1) %>%
    filter(CONSEQUENCE != "synonymous") %>% 
    filter(AD.NORMAL.1 != 0 | AD.NORMAL.2 != 0) %>%
    filter(AD.TUMOR.2 > 3) %>% 
    filter(FILTER == "PASS") %>% 
    filter(AF.NORMAL == 0) %>% 
    ungroup() %>% 
    mutate(GENEID = as.character(GENEID)) 
  
  org_translation <- select(org.Hs.eg.db, as.character(my_vcfs$GENEID), c("SYMBOL", "GENENAME"))
  
  my_vcfs <- left_join(my_vcfs, org_translation, by = c("GENEID" = "ENTREZID"))
  
  variants_w_o_pon <- my_vcfs %>%
    group_by(sample, snp_id) %>% 
    filter(row_number() == 1) %>%
    group_by(snp_id) %>% 
    summarise(recurrence = paste(as.character(sample), collapse = ";"))  %>% 
    merge(., my_vcfs, by = 'snp_id') %>% 
    group_by(sample, snp_id) %>% 
    filter(row_number() == 1) %>% 
    group_by(snp_id) %>%
    mutate(counts = n())  
     
  
  variants <- variants_w_o_pon %>% 
    filter(!grepl("N", recurrence)) %>% 
    ungroup() %>% 
    arrange(desc(counts))
  
  genes <- variants %>%
    group_by(GENEID) %>%
    mutate(gene_counts = sum(counts)) %>%
    mutate(gene_recurrence = paste(recurrence, collapse=";")) %>%
    mutate(gene_recurrence = map_chr(strsplit(as.character(gene_recurrence) ,";"), function(x) paste(unique(x), collapse=";"))) %>% 
    mutate(gene_recurrence_counts = map_int(strsplit(as.character(gene_recurrence) ,";"), function(x) length(unique(x)))) %>%
    dplyr::arrange(desc(counts)) %>% 
    ungroup() 
    
  
  samples <- variants %>% 
    ungroup() %>% 
    dplyr::select(sample, SYMBOL, GENENAME) %>% 
    dplyr::arrange(desc(sample)) 
    
    
  samples <- aggregate(SYMBOL ~ sample, data = unique(samples), paste, collapse = ",")
    
  my_list <-  list("variants_w_o_pon" = variants_w_o_pon, "variants" = variants, "genes" = genes, "samples" = samples)
    
  }

tidy_tumors <- retidy_vcfs(test0)
tidy_cell_lines <- retidy_vcfs(tidy_cell_lines0)
tidy_normals <- retidy_pon(tidy_normals0)
tidy_reynolds <- retidy_pon(tidy_reynolds0)

# find intersection between kooi and stachelek recurrent variants ---------

kooi_vars <-  read.table("doc/rb_variants_kooi_supp_info/tidy_format/variant_summary_kooi.csv") 
stchlk_kooi_var_intxn <-  inner_join(tidy_vcfs$variants, kooi_vars, by = c("seqnames" = "Chr", "start" = "Start", "end" = "End"))

# find intersection between kooi and stachelek recurrent genes ------------
kooi_genes <- read.table("doc/rb_variants_kooi_supp_info/tidy_format/gene_summary_kooi.csv")
stchlk_kooi_gene_intxn <- inner_join(tidy_vcfs$genes, kooi_genes, by = c("SYMBOL" = "gene")) %>% 
  arrange(desc(n_nonsilent))

tabulate_recurrence <-  function(recur_vec){
  table(unlist(strsplit(recur_vec, ";")))
}


# save tidy datasets as csvs ----------------------------------------------

save_tidy_csvs <- function(df_list, base_path){
  browser()
  dir.create(base_path)
  
  subpath_names <- names(df_list)
  
  subpaths <- lapply(subpath_names, function(x) paste0(base_path, x, "_tidy_table.csv"))

  clean_list_cols <- function(dataset2) {
    dataset2 <- dplyr::mutate(dataset2, recurrence = lapply(recurrence, paste, collapse = "; "), recurrence = unlist(recurrence))  
  }
  
  df_list[1:2] <- lapply(df_list[1:2], clean_list_cols)
  
  purrr::map2(df_list, subpaths, function(x,y) write.table(x, y, sep = ",", row.names = FALSE))
  
}

save_tidy_csvs(tidy_reynolds, "~/rb_pipeline/results/20170807_vcfs_reynolds/")


# read back in tidy csvs --------------------------------------------------

read_tidy_csvs <- function(basepath){
  my_path <- list.files(basepath, full.names = TRUE)
  tidy_vcfs <- lapply(my_path, read.table, sep=",", header  = TRUE)
  tidy_vcfs_names <- gsub("_.*$", "", basename(my_path))
  tidy_vcfs <- setNames(tidy_vcfs, tidy_vcfs_names)
  
}

tidy_reynolds <- read_tidy_csvs("~/rb_pipeline/results/20170807_vcfs_reynolds/")
tidy_cell_lines <- read_tidy_csvs("~/rb_pipeline/results/20170809_tidy_vcfs_cell_lines/") 
tidy_tumors<- read_tidy_csvs("~/rb_pipeline/results/20170809_tidy_vcfs_tumors/") 
