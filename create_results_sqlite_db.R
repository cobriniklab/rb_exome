library(DBI)
library(tidyverse)
library('rprojroot')
library('glue')
library('janitor')
library('zeallot')
library(fs)

proj_dir = rprojroot::find_root(criterion = has_file_pattern("*.Rproj"))

supplemental_db <- dbConnect(RSQLite::SQLite(), fs::path(proj_dir, "doc/RB_exome_manuscript/stachelek_supplemental/supplemental_tables.sqlite"))

supplemental_csvs <- fs::path(proj_dir, "doc/RB_exome_manuscript/stachelek_supplemental/") %>% 
  dir_ls(regex = ".*[0-9]+.*.csv") %>% 
  purrr::set_names(fs::path_file(.)) %>% 
  purrr::map(read_csv) %>% 
  identity()

scheme_table_names <- names(supplemental_csvs)

descriptive_table_names <- c("patients_in_study", "rb1_mycn_status_in_tumor_cell_line",
                             "rb_somatic_variants_all_studies", "sample_tally_from_ngs_studies",
                             "bcor_aberrations", "SCNA_in_rb_tumor_cell_line_by_gene",
                             "SCNA_in_rb_tumor_cell_line_by_sample", "coverage_by_study",
                             "somatic_variants_from_VC_tumor_cell_line_by_strelka_m2",
                             "somatic_variants_from_reynolds_cell_line_by_m2",
                             "str_typing_for_rb_cell_line", "allelic_imbalance_in_tumor_cell_line",
                             "variants_by_lawrence", "recurrent_somatic_variants_vc_tumor_cell_line",
                             "recurrent_somatic_variant_reynolds_cell_line", 
                             "days_in_culture_cell_line", "survey_vaf_at_recurrent_sites")

names(scheme_table_names) <- descriptive_table_names

safe_write <- safely(dbWriteTable)

my_res <- purrr::imap(supplemental_csvs, ~safe_write(supplemental_db, .y, .x))

dbListTables(supplemental_db)
  
dbDisconnect(supplemental_db)

# unlink("doc/RB_exome_manuscript/stachelek_supplemental/supplemental_tables.sqlite")


