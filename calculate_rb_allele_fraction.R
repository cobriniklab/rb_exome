# SCNA ------------------------------
rb_scna_hits <- 
  read_csv("stachelek_supplemental/table_s04.csv") %>% 
  dplyr::filter(SYMBOL == "RB1") %>%
  tidyr::pivot_longer(matches("^RB-*"), names_to = "sample_id", values_to = "SCNA") %>% 
  dplyr::filter(SCNA < 1.8) %>% 
  dplyr::select(sample_id, SCNA) %>%
  identity()

rb_scna_close_look <- read_csv("stachelek_supplemental/table_s05.csv") %>% 
  dplyr::filter(chrom == 13, sample_id == "CHLA-VC-RB33-T") %>% 
  # dplyr::filter(SYMBOL == "RB1") %>%
  # tidyr::pivot_longer(matches("^RB-*"), names_to = "sample_id", values_to = "SCNA") %>% 
  # dplyr::filter(SCNA < 1.8) %>% 
  # dplyr::select(sample_id, SCNA) %>%
  identity()


# SNV ------------------------------
rb_vaf_hits <- read_csv("stachelek_supplemental/table_02.csv") %>% 
  dplyr::filter() %>% 
  identity()

rb_minimal_allele_fraction <- tribble(
  ~Sample, ~`Minimal Allele Fraction`,
  "RB-14-T", 0.74,
  "RB-20-T", NA,
  "RB-24-T", 'MYCNA',
  "RB-28-T", 0.47,
  "RB-29-T", 0.99,
  "RB-31-T", 0.79,
  "RB-33-T", 0.87,
  "RB-41-T", NA,
  "RB-43-T", NA,
  "RB-46-T", 0.32,
  "RB-48-T", 0.47,
  "RB-49-T", 0.84
)

# total------------------------------
rb_hits <- dplyr::full_join(rb_vaf_hits, rb_scna_hits, by = c("Sample" = "sample_id")) %>% 
  dplyr::left_join(rb_minimal_allele_fraction, by = "Sample") %>%
  write_csv("stachelek_supplemental/table_02.csv") %>% 
  identity()
