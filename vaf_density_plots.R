library(patchwork)
library(cowplot)

# load data------------------------------
vc_somatic_variants_path <- fs::path("~/rb_pipeline/doc/RB_exome_manuscript/stachelek_supplemental", "table_s06.csv")
vc_somatic_variants <- read_csv(vc_somatic_variants_path)

reynolds_somatic_variants_path <- fs::path("~/rb_pipeline/doc/RB_exome_manuscript/stachelek_supplemental/table_s07.csv")
reynolds_somatic_variants <- read_csv(reynolds_somatic_variants_path)

# plots ------------------------------
m2_vars <- 
  vc_somatic_variants %>% 
  dplyr::filter(!is.na(AF.TUMOR.m2)) %>% 
  dplyr::mutate(VAF = AF.TUMOR.m2) %>% 
  dplyr::select(-AF.TUMOR.strelka, -AF.TUMOR.m2) %>% 
  dplyr::mutate(`variant caller` = "m2") %>% 
  identity()

strelka_vars <- 
  vc_somatic_variants %>% 
  dplyr::filter(!is.na(AF.TUMOR.strelka)) %>% 
  dplyr::mutate(VAF = AF.TUMOR.strelka) %>% 
  dplyr::select(-AF.TUMOR.strelka, -AF.TUMOR.m2) %>% 
  dplyr::mutate(`variant caller` = "strelka") %>% 
  identity()

all_vars <- dplyr::bind_rows(m2_vars, strelka_vars)  
  
vaf_dist_plot <- ggplot(all_vars, aes(x = VAF, color = `variant caller`)) + 
  geom_density() + 
  theme_cow_plot()

# ------------------------------

# tally variants per type

m2_vars %>% 
  dplyr::mutate(sample_type = str_extract(sample, "[A-Z]+$")) %>% 
  split(.$sample_type) %>% 
  identity() %>% 
  purrr::map(dim)

strelka_vars %>% 
  dplyr::mutate(sample_type = str_extract(sample, "[A-Z]+$")) %>% 
  split(.$sample_type) %>% 
  identity() %>% 
  purrr::map(dim) %>% 
  identity()

# check concordance------------------------------

m2_vars %>% 
  dplyr::filter(VAF > 0.05) %>% 
  dim()

strelka_vars %>% 
  dplyr::filter(VAF > 0.05) %>% 
  dim()

strelka_v_m2 <- anti_join(strelka_vars, m2_vars, by = c("sample", "SYMBOL"))

m2_v_strelka <- anti_join(m2_vars, strelka_vars, by = c("sample", "SYMBOL"))

# venn overlap ------------------------------
library(eulerr)

# format tibble
test0 <- 
  vc_somatic_variants %>% 
  dplyr::select(SYMBOL, VAF.m2 = AF.TUMOR.m2, VAF.strelka = AF.TUMOR.strelka) %>% 
  dplyr::mutate(VAF.strelka = as.logical(VAF.strelka)) %>% 
  dplyr::mutate(VAF.m2 = as.logical(VAF.m2)) %>% 
  identity()

test0[is.na(test0)] <- FALSE

fit <- euler(test0[, -1])

venn_caller_overlap <- plot(fit)

out_plot <- vaf_dist_plot + 
  venn_caller_overlap + 
  plot_layout(widths = c(2, 1)) + 
  plot_annotation(tag_levels = 'A')

ggsave("~/rb_pipeline/doc/RB_exome_manuscript/stachelek_supplemental/fig_s02.pdf", out_plot)
