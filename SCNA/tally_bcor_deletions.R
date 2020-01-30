vaf_cols <- c("VAF.m2", "VAF.strelka")

kooi_bcor_deletions = c("T22", "T67")

zhang_bcor_deletions <- 

test0 <- 
  recurrent_reynolds_cl_vars %>% 
  # dplyr::filter(!is.na(recurrent_sample)) %>% 
  dplyr::mutate_at(vars(vaf_cols), as.numeric) %>% 
  dplyr::mutate(VAF = coalesce(VAF.m2, VAF.strelka)) %>% 
  dplyr::filter(VAF.m2 > 0.45) %>%
  identity()
  
  