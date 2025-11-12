# this script contain the function to extract the latest assay status for a uniprot ID
library(writexl)

get_assay_status <- function(uniprot_list, 
                             assay_categories,
                             antigen_categories,
                             in_product,
                             d72,
                             master_list){
  
  assay_status_categories_missing <-
    d72 %>%
    select(assay_status) %>%
    filter(!is.na(assay_status)) %>%
    mutate(assay_status = str_to_lower(assay_status)) %>%
    distinct() %>%
    filter(!assay_status %in% assay_categories$assay_status)
  
  antigen_status_categories_missing <-
    master_list %>%
    pivot_longer(cols = c(absea, icosagen, bon_opus),
                 names_to = "vendor_antigen",
                 values_to = "vendor_categry") %>%
    mutate(vendor_categry = tolower(vendor_categry)) %>%
    select( vendor_categry) %>%
    distinct() %>%
    filter(!is.na(vendor_categry)) %>%
    filter(!vendor_categry %in% antigen_categories$antigen_status)
  
  if(assay_status_categories_missing %>% nrow() > 0 |
     antigen_status_categories_missing %>% nrow() > 0){
    
    cat("⚠️ Missing categories detected! ⚠️\n")
    
    if (nrow(assay_status_categories_missing) > 0) {
      cat("\n🔍 Missing assay status categories:\n")
      print(assay_status_categories_missing %>% pull(assay_status))
    }
    
    if (nrow(antigen_status_categories_missing) > 0) {
      cat("\n🔍 Missing antigen status categories:\n")
      print(antigen_status_categories_missing %>% pull(vendor_categry))
    }
    
    cat("\n📋 Please update the 'status_categories_assay_antigen.xlsx' file with the missing categories above.\n")
    
  } else{
    
    status_level_adjust <-
      uniprot_list %>%
      left_join(d72) %>%
      mutate(assay_status = if_else(is.na(assay_status), "not screened", assay_status)) %>%
      mutate(assay_status = str_to_lower(assay_status)) %>%
      left_join(assay_categories) %>%
      mutate(assay_status_category = if_else(assay_status == "not screened", "Not screened", assay_status_category),
             antigen_check_needed = if_else(assay_status == "not screened", "yes", antigen_check_needed)) %>%
      mutate(assay_order_in_pipeline_same_ab = case_when(assay_status_category == "Terminated" ~ 15,
                                                         .default = assay_order_in_pipeline)) %>%
      arrange(desc(assay_order_in_pipeline_same_ab)) %>%
      group_by(uniprot, vendor_a) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      arrange(desc(assay_order_in_pipeline)) %>%
      mutate(assay_status = case_when((!is.na(assay_synonym_to)) ~ assay_synonym_to,
                                      .default = assay_status)) %>%
      left_join(in_product) %>%
      mutate(assay_status_category = case_when(!(is.na(in_product)) ~ in_product,
                                               .default = assay_status_category)) %>%
      select(-in_product, -assay_synonym_to)
    
    status_level_1 <-
      status_level_adjust %>%
      group_by(uniprot) %>%
      slice_head(n = 1) %>%
      ungroup()
    
    status_screening_clear <-
      status_level_1 %>%
      filter(antigen_check_needed == "no" | assay_status_category %in% c("In Maximus", "In released product")) %>%
      select(uniprot, assay_status_category, assay_status, art_no_ag, artnr_a_arm, vendor_a, art_nr_b_arm, vendor_b, antigen, vendor_ag, issue, comment_status)
    
    status_antigen <-
      status_level_adjust %>%
      filter(!uniprot %in% status_screening_clear$uniprot) %>%
      left_join(master_list %>%
                  select(uniprot, absea, icosagen, "bonopus" = bon_opus), by = "uniprot") %>%
      select(names(status_level_1), absea, icosagen, bonopus, art_no_ag) %>%
      mutate(across(c(vendor_a, absea, icosagen, bonopus, art_no_ag), ~ str_to_lower(.x))) %>%
      mutate(across(c(vendor_a, absea, icosagen, bonopus), ~ str_replace_all(.x, "0", NA_character_)))
    
    never_initiated <-
      status_antigen %>%
      filter(assay_status_category == "Not screened" &
               is.na(absea) &
               is.na(icosagen) &
               is.na(bonopus)) %>%
      mutate(antigen_status_category = "Never initiated",
             antigen_status = "Never initiated") %>%
      select(uniprot, assay_status_category, assay_status, antigen_status_category, art_no_ag, artnr_a_arm, vendor_a, art_nr_b_arm, vendor_b, antigen, vendor_ag, issue, comment_status)
    
    not_screened_ag_initiated <-
      status_antigen %>%
      filter(assay_status_category == "Not screened" & !uniprot %in% never_initiated$uniprot) %>%
      pivot_longer(cols = c(absea, icosagen, bonopus),
                   names_to = "vendor_antigen",
                   values_to = "vendor_categry") %>%
      left_join(antigen_categories, by = c("vendor_categry" = "antigen_status")) %>%
      filter(!is.na(vendor_categry)) %>%
      arrange(desc(antigen_order_in_pipeline)) %>%
      group_by(uniprot) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      select(names(never_initiated))
    
    screened_agrisera_prep <-
      status_antigen %>%
      filter(str_detect(vendor_a, "agrisera")) %>% 
      mutate(vendor_a = str_remove(vendor_a, "/abcore")) %>%
      mutate(vendor_a = case_when(vendor_a == "agrisera" ~ paste0("agrisera/", art_no_ag),
                                  .default = vendor_a)) %>%
      mutate(vendor_a = str_remove_all(vendor_a, "agrisera")) %>%
      mutate(vendor_a = str_remove_all(vendor_a, "/")) %>%
      pivot_longer(cols = c(absea, icosagen, bonopus),
                   names_to = "vendor_antigen",
                   values_to = "vendor_categry") %>%
      left_join(antigen_categories, by = c("vendor_categry" = "antigen_status"))
    
    screened_agrisera_1 <-
      screened_agrisera_prep %>%
      filter(!is.na(vendor_categry)) %>%
      filter(vendor_a != vendor_antigen) %>%
      arrange(desc(antigen_order_in_pipeline)) %>%
      group_by(uniprot, vendor_a) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      select(names(never_initiated))
    
    screened_agrisera_2 <-
      screened_agrisera_prep %>%
      filter(is.na(vendor_categry)) %>%
      filter(!uniprot %in% screened_agrisera_1$uniprot) %>%
      arrange(desc(assay_order_in_pipeline)) %>%
      group_by(uniprot) %>%
      slice_head(n = 1) %>%
      mutate(antigen_status_category = "Never re-initiated") %>%
      select(names(never_initiated))
    
    screened_agrisera <-
      bind_rows(screened_agrisera_1,
                screened_agrisera_2)
      
    screened_other <-
      status_antigen %>%
      filter(!uniprot %in% never_initiated$uniprot &
               !uniprot %in% not_screened_ag_initiated$uniprot &
               !uniprot %in% screened_agrisera$uniprot) %>%
      pivot_longer(cols = c(absea, icosagen, bonopus),
                   names_to = "vendor_antigen",
                   values_to = "vendor_categry") %>%
      left_join(antigen_categories, by = c("vendor_categry" = "antigen_status")) %>%
      arrange(desc(antigen_order_in_pipeline)) %>%
      group_by(uniprot) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      mutate(antigen_status_category = if_else(is.na(antigen_status_category), "Never initiated", antigen_status_category)) %>%
      select(names(never_initiated))
    
    status_antigen_final <-
      bind_rows(never_initiated,
                not_screened_ag_initiated,
                screened_agrisera,
                screened_other)
    
    status_final <-
      bind_rows(status_screening_clear %>%
                  mutate(antigen_status_category = NA),
                status_antigen_final) %>%
      mutate(assay_status_antigen_status = case_when(is.na(antigen_status_category) ~ assay_status_category,
                                                     .default = paste0(assay_status_category, ";", antigen_status_category))) %>%
      relocate(antigen_status_category, .after = assay_status) %>%
      relocate(assay_status_antigen_status, .after = uniprot)
    
    write_xlsx(status_final, paste0("output/", Sys.Date(), "_assay_status_per_uniprot.xlsx"))
    
    return(status_final)
  }
}