
library(stchlk.exome)
library(dplyr)

ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")

#clinvar <-  read.table('../Homo_sapiens/ClinVar/clinvar_alleles.single.b37.tsv',sep='\t',header=T,quote='',comment.char='')


# load input --------------------------------------------------------------

mutect2_tumor_filenames <- list.files(path="./output/gatk", pattern="[[:digit:]]-T_1_.*anno.vcf.gz$", full.names = TRUE)
mutect2_cell_line_filenames <- list.files(path="./output/gatk", pattern="[[:digit:]]-CL_1_.*anno.vcf.gz$", full.names = TRUE)
mutect2_pon_filenames <- list.files(path="./output/gatk", pattern="[[:digit:]]-N_1_.*anno.vcf.gz$", full.names = TRUE)

strelka_snv_filenames <- list.files(path="./output/strelka/", pattern=".*somatic.snvs.vcf.gz$",  full.names = TRUE, recursive = TRUE)
strelka_indel_filenames <- list.files(path="./output/annovar/", pattern="somatic.indels.vcf.gz$",  full.names = TRUE, recursive = TRUE)

rb_filenames <- list.files(path="./output/gatk", pattern="*rb_anno.vcf.gz$", full.names = TRUE)
pon_rb_filenames <-  list.files(path="./output/annovar", pattern="*pon_rb_annovar.vcf.gz$", full.names = TRUE)

mutect2_tumor_names <-substr(mutect2_tumor_filenames,15,18)
mutect2_cell_line_names <-substr(mutect2_cell_line_filenames,15,18)
mutect2_normal_names <-substr(mutect2_pon_filenames,15,18)

strelka_snv_names <- substr(strelka_snv_filenames ,19,22)
strelka_indel_names <- substr(strelka_indel_filenames ,19,22)

rb_names <- substr(rb_filenames, 15, 18)
pon_rb_names <- substr(pon_rb_filenames, 18, 21)

# read vcf ----------------------------------------------------------------


#svp <- ScanVcfParam(info="AF", geno="GT")

mutect2_tumor_list <- lapply(mutect2_tumor_filenames, function(x)readVcf(x, "hg19"))
mutect2_cell_line_list <- lapply(mutect2_cell_line_filenames, function(x)readVcf(x, "hg19"))
mutect2_pon_list <- lapply(mutect2_pon_filenames, function(x)readVcf(x, "hg19"))

strelka_snv_list <- lapply(strelka_snv_filenames, read.table, sep = ",", header = TRUE)
strelka_indel_list <- lapply(strelka_indel_filenames, function(x)readVcf(x, "hg19"))

pon_rb_list <- map(pon_rb_filenames, function(x)readVcf(x,"hg19"))

fl <- mutect2_tumor_filenames[[1]]

names(mutect2_tumor_list) <- mutect2_tumor_names
names(mutect2_cell_line_list) <- mutect2_cell_line_names
names(mutect2_pon_list) <- mutect2_normal_names

names(strelka_snv_list) <- strelka_snv_names
names(strelka_indel_list) <- strelka_indel_names

saveRDS(strelka_snv_list, "doc/strelka_snv_list.rds")

strelka_snv_list <- readRDS("doc/strelka_snv_list.rds")





strelka_snvs <- map_dfr(strelka_snv_list, clean_strelka_csv, .id = "Sample")
saveRDS(strelka_snvs, "doc/strelka_snvs.rds")
strelka_snvs <- readRDS("doc/strelka_snvs.rds")




strelka_snvs0 <- strelka_snvs %>%
  rowwise() %>%
  mutate(TUMOR_VAF = calc_strelka_VAF(c(`TUMOR-AU`,`TUMOR-CU`,`TUMOR-GU`,`TUMOR-TU`))) %>%
  mutate(NORMAL_VAF = calc_strelka_VAF(c(`NORMAL-AU`,`NORMAL-CU`,`NORMAL-GU`,`NORMAL-TU`)))

saveRDS(strelka_snvs0, "doc/strelka_snvs0.rds")


names(pon_rb_list) <- pon_rb_names

keep_names <- c("24-T", "28-C", "28-T", "29-T", "29-C", "31-T", "31-C", "41-C", "46-T", "46-C", "48-T", "48-C")
rb_list <- rb_list[names(rb_list) %in% keep_names] #the set of rb-specific vcfs that have any variants

vcf <-  mutect2_tumor_list[[1]]

test_tumor_list <- lapply(mutect2_tumor_list, '[', c(1:10))
test_reynolds_list <- lapply(mutect2_reynolds_list, '[', c(1:10))

tidy_m2_tumors0 <- collate_vcfs(mutect2_tumor_list, mutect2_tidy)
tidy_m2_cell_lines0 <- collate_vcfs(mutect2_cell_line_list, mutect2_tidy)
tidy_m2_normals0 <- collate_pon(mutect2_pon_list, mutect2_tidy)
tidy_m2_rb0 <- collate_vcfs(rb_list, mutect2_tidy)
tidy_m2_pon_rb0 <- collate_pon_rbs(pon_rb_list, mutect2_tidy)

tidy_strelka_snvs0 <- collate_vcfs(strelka_snv_list, strelka_tidy)



# saveRDS(tidy_rb0, "./results/SNV/tidy_mutect2_rb1_variants.rds")
saveRDS(tidy_tumors0, "./results/SNV/tidy_mutect2_anno_tumor_variants_2nd_ed.rds")
saveRDS(tidy_normals0, "./results/SNV/tidy_mutect2_anno_pon_variants_2nd_ed.rds")
saveRDS(tidy_cell_lines0, "./results/SNV/tidy_mutect2_anno_cell_line_variants_2nd_ed.rds")
saveRDS(tidy_rb0, "./results/SNV/tidy_mutect2_anno_rb_variants_2nd_ed.rds")
saveRDS(tidy_strelka_snvs0, "./results/SNV/tidy_strelka_snvs.rds")
# saveRDS(tidy_pon_rb0, "./results/SNV/tidy_mutect2_anno_pon_rb_variants.rds")

tidy_tumors0 <- readRDS("./results/SNV/tidy_mutect2_anno_tumor_variants_2nd_ed.rds")
tidy_cell_lines0 <- readRDS("./results/SNV/tidy_mutect2_anno_cell_line_variants_2nd_ed.rds")
tidy_normals0 <- readRDS("./results/SNV/tidy_mutect2_anno_pon_variants_2nd_ed.rds")


# retidy vcfs -------------------------------------------------------------

# retidy_pon <- function(my_pon){
# # germline_vcfs <- function(my_pon){
#   browser()
#   my_pon <- dplyr::arrange(my_pon, desc(AF.NORMAL)) %>%
#     dplyr::group_by(sample, snp_id) %>%
#     dplyr::filter(row_number() == 1) %>%
#     dplyr::group_by(snp_id) %>%
#     dplyr::mutate(recurrence=paste(sample,collapse=';')) %>%
#     dplyr::mutate(counts = n()) %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(Gene.refGene = as.character(Gene.refGene))
#
#   variants <- my_pon %>%
#     group_by(snp_id) %>%
#     filter(AF.NORMAL > 0.4) %>%
#     filter(row_number() == 1) %>%
#     ungroup() %>%
#     arrange(desc(counts))
#
#   genes <- variants %>%
#     group_by(snp_id) %>%
#     filter(row_number() == 1) %>%
#     group_by(Gene.refGene) %>%
#     mutate(gene_counts = sum(counts)) %>%
#     mutate(gene_recurrence = paste(recurrence, collapse=";")) %>%
#     mutate(gene_recurrence = map_chr(strsplit(as.character(gene_recurrence) ,";"), function(x) paste(unique(x), collapse=";"))) %>%
#     mutate(gene_recurrence_counts = map_int(strsplit(as.character(gene_recurrence) ,";"), function(x) length(unique(x)))) %>%
#     dplyr::arrange(desc(counts)) %>%
#     ungroup()
#
#
#   samples <- variants %>%
#     ungroup() %>%
#     dplyr::select(sample, Gene.refGene) %>%
#     dplyr::arrange(desc(sample))
#
#
#   samples <- aggregate(Gene.refGene ~ sample, data = unique(samples), paste, collapse = ",")
#
#   my_list <-  list("variants" = variants, "genes" = genes, "samples" = samples)
#
#   }


tidy_tumors <- retidy_vcfs(tidy_tumors0, tidy_normals0)
tidy_cell_lines <- retidy_vcfs(tidy_cell_lines0, tidy_normals0)

# tidy_normals <- germline_vcfs(tidy_normals0)

# tidy_pon_rb <-  group_by(tidy_normals_tst, seqnames, start, end) %>%
#   mutate(recurrence = paste(sample, collapse = ";")) %>%
#   mutate(recurrence = map_chr(strsplit(as.character(recurrence) ,";"), function(x) paste(unique(x), collapse=";"))) %>%
#   mutate(recurrence_counts = map_int(strsplit(as.character(recurrence) ,";"), function(x) length(unique(x)))) %>%
#   arrange(recurrence_counts) %>%
#   filter(FILTER == "PASS")
#
# write.table(tidy_pon_rb, "./results/SNV/tidy_pon_rb.csv", sep = ",", row.names = FALSE)


stchlk.exome::save_tidy_csvs(tidy_reynolds, "~/rb_pipeline/results/20170807_tidyvcfs_reynolds/")
stchlk.exome::save_tidy_csvs(tidy_cell_lines, "~/rb_pipeline/results/SNV/20170921_tidy_vcfs_cell_lines/")
stchlk.exome::save_tidy_csvs(tidy_tumors, "~/rb_pipeline/results/SNV/20170921_tidy_vcfs_cell_tumors/")

tidy_reynolds <- stchlk.exome::read_tidy_csvs("~/rb_pipeline/results/SNV/20170921_tidy_vcfs_reynolds")
tidy_cell_lines <- stchlk.exome::read_tidy_csvs("~/rb_pipeline/results/SNV/20170921_tidy_vcfs_cell_lines")
tidy_tumors<- stchlk.exome::read_tidy_csvs("~/rb_pipeline/results/SNV/20170921_tidy_vcfs_tumors")




# find intersection between tumor and cell line variants ------------------
t_cl_vars <-  inner_join(tidy_tumors$variants, tidy_cell_lines$variants, by = c("seqnames", "start", "end", "REF", "ALT")) %>%
  mutate(t_prefix = substr(sample.x, 1,2)) %>%
  mutate(cl_prefix = substr(sample.y, 1,2)) %>%
  filter(t_prefix == cl_prefix | is.na(t_prefix) | is.na(cl_prefix))%>%
  filter(AF.TUMOR.x > 0.279 | AF.TUMOR.y > 0.279) %>%
  arrange(Gene.refGene.x, sample.x) %>%
  ungroup() %>%
  dplyr::select(seqnames, start, Gene.refGene = Gene.refGene.x, sample.Tumor = sample.x,
                AF.Tumor = AF.TUMOR.x, sample.Cell.Line = sample.y, AF.Cell.Line = AF.TUMOR.y,
                Func.refGene.value = Func.refGene.value.y, GeneDetail.refGene = GeneDetail.refGene.y,
                ExonicFunc.refGene = ExonicFunc.refGene.y, AAChange.refGene = AAChange.refGene.y) %>%
  mutate(AF.Tumor = as.numeric(AF.Tumor), AF.Cell.Line = as.numeric(AF.Cell.Line))



t_cl_vars_besides_RB1 <- filter(t_cl_vars, Gene.refGene != "RB1")

write.table(t_cl_vars, "./results/SNV/t_cl_vars.csv", sep =",", row.names = FALSE)
write.table(t_cl_vars_besides_RB1, "./results/SNV/t_cl_vars_besides_RB1.csv", sep =",", row.names = FALSE)

# find intersection between zhang and stachelek recurrent genes -----------

stchlk_overall_variants <- dplyr::select(tidy_tumors$variants, -33, -34, -35, -36) %>%
  rbind(tidy_cell_lines$variants) %>%
  ungroup()

zhang_genes <- read.table("doc/zhang_supp_info/tidy_format/variants_summary_pdf.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

zhang_genes <- dplyr::rename(zhang_genes, Gene.refGene = Gene) %>%
  group_by(Gene.refGene) %>%
  mutate(zhang_incidence = n()) %>%
  filter(row_number() == 1)

stchlk_zhang_gene_intxn <-  inner_join(stchlk_overall_variants, zhang_genes, by = "Gene.refGene") %>%
  arrange(Gene.refGene, sample) %>%
  mutate(sample.Tumor = ifelse(grepl("T", sample), as.character(sample), "")) %>%
  mutate(sample.Cell.Line = ifelse(grepl("C", sample), as.character(sample),"")) %>%
  mutate(AF.Tumor = ifelse(grepl("T", sample), as.character(AF.TUMOR), "")) %>%
  mutate(AF.Cell.Line = ifelse(grepl("C", sample), as.character(AF.TUMOR), "")) %>%
  dplyr::select(Gene.refGene, sample.Tumor, AF.Tumor, sample.Cell.Line, AF.Cell.Line,
                Func.refGene.value, GeneDetail.refGene, ExonicFunc.refGene, AAChange.refGene, zhang_incidence) %>%
  filter(Gene.refGene != "RB1") %>%
  mutate(AF.Tumor = as.numeric(AF.Tumor), AF.Cell.Line = as.numeric(AF.Cell.Line))

write.table(stchlk_zhang_gene_intxn, "./doc/SNV/stchlk_zhang_gene_intxn.csv", sep =",", row.names = FALSE)
# only BCOR!

# find intersection between kooi and stachelek recurrent variants ---------

kooi_vars <-  read.table("doc/rb_variants_kooi_supp_info/tidy_format/variant_summary_kooi.csv")
kooi_vars <- dplyr::rename(kooi_vars, seqnames = Chr, start = Start, end = End, REF = Ref, ALT = Obs) %>%
  group_by(Gene) %>%
  mutate(kooi_incidence = n()) %>%
  filter(row_number() == 1)
stchlk_kooi_var_intxn <- semi_join(tidy_vcfs$variants, kooi_vars, by = c("seqnames", "start", "end", "REF", "ALT"))
write.table(stchlk_kooi_var_intxn, "./doc/SNV/stchlk_kooi_var_intxn.csv", sep =",", row.names = FALSE)

# find intersection between kooi and stachelek recurrent genes ------------
stchlk_kooi_gene_intxn <-  inner_join(stchlk_overall_variants, kooi_vars, by = c("Gene.refGene" = "Gene")) %>%
  arrange(Gene.refGene, sample) %>%
  mutate(sample.Tumor = ifelse(grepl("T", sample), as.character(sample), "")) %>%
  mutate(sample.Cell.Line = ifelse(grepl("C", sample), as.character(sample),"")) %>%
  mutate(AF.Tumor = ifelse(grepl("T", sample), AF.TUMOR, "")) %>%
  mutate(AF.Cell.Line = ifelse(grepl("C", sample), AF.TUMOR, ""))  %>%
  dplyr::select(seqnames = seqnames.x, start = start.x, Gene.refGene, sample.Tumor, AF.Tumor, sample.Cell.Line, AF.Cell.Line,
                Func.refGene.value, GeneDetail.refGene, ExonicFunc.refGene, AAChange.refGene, kooi_incidence) %>%
  filter(Gene.refGene != "RB1") %>%
  mutate(AF.Tumor = as.numeric(AF.Tumor), AF.Cell.Line = as.numeric(AF.Cell.Line))


write.table(stchlk_kooi_gene_intxn, "./doc/SNV/stchlk_kooi_gene_intxn.csv", sep =",", row.names = FALSE)

t_cl_vars_mutect2 <- t_cl_vars_besides_RB1 %>%
  dplyr::bind_rows(stchlk_zhang_gene_intxn) %>%
  dplyr::bind_rows(stchlk_kooi_gene_intxn) %>%
  group_by(Gene.refGene, sample.Tumor) %>%
  mutate(seqnames = as.character(seqnames)) %>%
  summarise_all(max, na.rm=TRUE) %>%
  group_by(Gene.refGene, sample.Cell.Line) %>%
  summarise_all(max, na.rm=TRUE)

write.table(t_cl_vars_mutect2, "./results/SNV/mutect2_t_cl_vars.csv", sep =",", row.names = FALSE)
t_cl_vars_mutect2 <- read.table("results/SNV/mutect2_t_cl_vars.csv", sep = ",", header = TRUE)

# compare mutect2 with strelka calls --------------------------------------

test <- inner_join(tidy_tumors$variants, strelka_variants$strelka_snvs, by = c("start" = "POS", "sample" = "Sample")) %>%
  arrange(Gene.refGene)

test0 <- inner_join(t_cl_vars_mutect2, strelka_variants$strelka_snvs, by = c("start" = "POS", "sample.Tumor" = "Sample", "sample.Cell.Line" = "Sample")) %>%
  arrange(Gene.refGene) %>%
  dplyr::select(-c(start.y:LRT_pred)) %>%
  mutate_if(is.factor, as.character) %>%
  group_by(start, sample.Tumor)

test2 <- inner_join(t_cl_vars_mutect2, strelka_variants$strelka_indels, by = c("start" = "POS", "sample.Tumor" = "Sample", "sample.Cell.Line" = "Sample")) %>%
  arrange(Gene.refGene) %>%
  dplyr::select(-c(start.y:LRT_pred)) %>%
  mutate_if(is.factor, as.character) %>%
  group_by(start, sample.Tumor)

stchlk_overall_variants0 <- dplyr::select(stchlk_overall_variants, c(sample, seqnames, start, Gene.refGene))

test3 <- inner_join(stchlk_overall_variants0, strelka_variants$strelka_indels, by = c("start" = "POS", "sample" = "Sample")) %>%
  arrange(Gene.refGene) %>%
  dplyr::select(-c(start.y:LRT_pred)) %>%
  mutate_if(is.factor, as.character)

test4 <- inner_join(stchlk_overall_variants0, strelka_variants$strelka_snvs, by = c("start" = "POS", "sample" = "Sample")) %>%
  arrange(Gene.refGene) %>%
  dplyr::select(-c(start.y:LRT_pred)) %>%
  mutate_if(is.factor, as.character)

unf_t_cl_vars <- rbind(test3, test4) %>%
  arrange(Gene.refGene)

t_cl_vars_manuscript <- rbind(test0, test2) %>%
  filter(Func.refGene.value != "intronic") %>%
  arrange(Gene.refGene)

write.table(t_cl_vars_manuscript, "./doc/SNV/t_cl_vars.csv", sep =",", row.names = FALSE)
t_cl_vars_manuscript <- read.table("doc/SNV/t_cl_vars.csv", sep = ",", header = TRUE)


# play with somatic signatures --------------------------------------------

read_vr <- function(vr0, sample_name){
  vr <- as(test, "VRanges")
  
  #Only SNV substitutions are currently supported; need to remove indels
  idx_snv = ref(vr) %in% DNA_BASES & alt(vr) %in% DNA_BASES
  vr = vr[idx_snv]
  vr$sample_name <- sample_name
  
  #drop softFilterMatrix, which would otherwise summarize GATK filter field, disallow concat of mult. vranges
  resetFilter(vr)
  return(vr)
}

vr_list <- map2(mutect2_tumor_list, names(mutect2_tumor_list), read_vr)

vr0 <- Reduce(c, vr_list)

sca_motifs = mutationContext(vr0, BSgenome.Hsapiens.UCSC.hg19)

sca_mm = motifMatrix(sca_motifs, group = "sample_name", normalize = TRUE)

head(round(sca_mm, 4))

pdf("./results/SNV/somatic_signatures_stachelek.pdf")
plotMutationSpectrum(sca_motifs, "sample_name")
dev.off()

n_sigs = 2:8

gof_nmf = assessNumberSignatures(sca_mm, n_sigs, nReplicates = 5)

gof_pca = assessNumberSignatures(sca_mm, n_sigs, pcaDecomposition)

plotNumberSignatures(gof_nmf)

plotNumberSignatures(gof_pca)

head(sca_motifs)




# plot VAF distribution  --------------------------------------------------

ggplot(tidy_tumors$variants, aes(AF.TUMOR)) + geom_histogram(binwidth = 0.01) +
  scale_x_continuous(breaks = pretty(tidy_tumors$variants$AF.TUMOR, n = 20)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(text = element_text(size=20))


ggplot(kooi_vars, aes(VAF)) + geom_histogram(binwidth = 0.01) +
  scale_x_continuous(breaks = pretty(tidy_tumors$variants$AF.TUMOR, n = 20)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(text = element_text(size=20))
