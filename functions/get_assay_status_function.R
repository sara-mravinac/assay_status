# this script contain the function to extract the latest assay status for a uniprot ID
library(writexl)

get_assay_status <- function(uniprot_list, 
                             assay_categories,
                             antigen_categories,
                             in_product,
                             d72,
                             master_list,
                             homology_list){
  
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
                 values_to = "vendor_category") %>%
    mutate(vendor_category = tolower(vendor_category)) %>%
    select( vendor_category) %>%
    distinct() %>%
    filter(!is.na(vendor_category)) %>%
    filter(!vendor_category %in% antigen_categories$antigen_status)
  
  if(assay_status_categories_missing %>% nrow() > 0 |
     antigen_status_categories_missing %>% nrow() > 0){
    
    cat("⚠️ Missing categories detected! ⚠️\n")
    
    if (nrow(assay_status_categories_missing) > 0) {
      cat("\n🔍 Missing assay status categories:\n")
      print(assay_status_categories_missing %>% pull(assay_status))
    }
    
    if (nrow(antigen_status_categories_missing) > 0) {
      cat("\n🔍 Missing antigen status categories:\n")
      print(antigen_status_categories_missing %>% pull(vendor_category))
    }
    
    cat("\n📋 Please update the 'status_categories_assay_antigen.xlsx' file with the missing categories above.\n")
    
  } else{
    
    if(d72 %>% filter(is.na(assay_status)) %>% nrow() > 0 ) cat("⚠️ The d72 file contains rows with missing 'status'! This can affect the results even if the script runs successfully. ⚠️\n\n")
    
    suppressMessages({
      
      # start from uniprot list and join with d72 assay data
      # add "not screened" status to uniprots that do not have data in d72
      # for the same assay (same uniprot, same vendor_a, same vendor_b, same artnr_a_arm, same art_nr_b_arm), keep only row with latest status
      # add information about assays in products, if status category disagree, in product overwrites any other status
      status_level_adjust <-
        uniprot_list %>%
        left_join(d72 %>%
                    mutate(assay_status = if_else(is.na(assay_status), "Unclear", assay_status))) %>%
        mutate(assay_status = if_else(is.na(assay_status), "not screened", assay_status)) %>%
        mutate(assay_status = str_to_lower(assay_status)) %>%
        left_join(assay_categories) %>%
        mutate(assay_status_category = case_when(assay_status == "not screened" ~ "Not screened",
                                                 assay_status == "unclear" ~ "Unclear", 
                                                 .default = assay_status_category),
               antigen_check_needed = if_else(assay_status == "not screened", "yes", antigen_check_needed)) %>%
        mutate(assay_order_in_pipeline_same_ab = case_when(assay_status_category == "Terminated" ~ 15,
                                                           .default = assay_order_in_pipeline)) %>%
        arrange(desc(assay_order_in_pipeline_same_ab)) %>%
        group_by(uniprot, vendor_a, vendor_b, artnr_a_arm, art_nr_b_arm) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        arrange(desc(assay_order_in_pipeline)) %>%
        mutate(assay_status = case_when((!is.na(assay_synonym_to)) ~ assay_synonym_to,
                                        .default = assay_status)) %>%
        left_join(in_product) %>%
        mutate(assay_status_category = case_when(!(is.na(in_product)) ~ in_product,
                                                 .default = assay_status_category)) %>%
        select(-in_product, -assay_synonym_to)
      
      # Separate those with clear screening status from those needing antigen status evaluation - assays with the highest order in pipeline is selected for the evaluation
      status_screening_clear <-
        status_level_adjust %>%
        arrange(desc(assay_order_in_pipeline)) %>%
        group_by(uniprot) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        filter(antigen_check_needed == "no" | assay_status_category %in% c("In Maximus", "In released product")) %>%
        select(uniprot, assay_status_category, assay_status, everything())
      
      # For the rest, add antigen data, all different assays for the same uniprot are kept
      status_antigen <-
        status_level_adjust %>%
        filter(!uniprot %in% status_screening_clear$uniprot) %>%
        left_join(master_list %>%
                    select(uniprot, absea, icosagen, "bonopus" = bon_opus), by = "uniprot") %>%
        select(names(status_level_adjust), absea, icosagen, bonopus, art_no_ag) %>%
        mutate(across(c(vendor_a, absea, icosagen, bonopus, art_no_ag), ~ str_to_lower(.x))) %>%
        mutate(across(c(vendor_a, absea, icosagen, bonopus), ~ str_replace_all(.x, "0", NA_character_)))
      
      # separate those never initiated on both screening and antigen level
      never_initiated <-
        status_antigen %>%
        filter(assay_status_category == "Not screened" &
                 is.na(absea) &
                 is.na(icosagen) &
                 is.na(bonopus)) %>%
        mutate(antigen_status_category = "Never initiated",
               antigen_status = "Never initiated") %>%
        select(uniprot, assay_status_category, assay_status, antigen_status_category, art_no_ag, artnr_a_arm, vendor_a, art_nr_b_arm, vendor_b, antigen, vendor_ag, issue, comment_status, sorter)
      
      # separate those never screened but with initiated antigen
      not_screened_ag_initiated <-
        status_antigen %>%
        filter(assay_status_category == "Not screened" & !uniprot %in% never_initiated$uniprot) %>%
        pivot_longer(cols = c(absea, icosagen, bonopus),
                     names_to = "vendor_antigen",
                     values_to = "vendor_category") %>%
        left_join(antigen_categories, by = c("vendor_category" = "antigen_status")) %>%
        filter(!is.na(vendor_category)) %>%
        arrange(desc(antigen_order_in_pipeline)) %>%
        group_by(uniprot) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        select(names(never_initiated))
      
      # screened which need antigen check
      screened <-
        status_antigen %>%
        filter(!uniprot %in% never_initiated$uniprot &
                 !uniprot %in% not_screened_ag_initiated$uniprot)
      
      # make sure that Agrisera is always in vendor_a position in case one of the vendors is agrisera and the other is not
      screened_switched_vendor <-
        screened %>%
        filter(str_detect(vendor_b, ("agrisera|Agrisera")) & !str_detect(vendor_a, ("agrisera|Agrisera"))) %>%
        mutate(vendor_b_temp = vendor_b,
               art_nr_b_arm_temp = art_nr_b_arm) %>%
        mutate(vendor_b = vendor_a,
               vendor_a = vendor_b_temp,
               art_nr_b_arm = artnr_a_arm,
               artnr_a_arm = art_nr_b_arm_temp) %>%
        select(-vendor_b_temp, -art_nr_b_arm_temp)
      
      # update screened with switched vendor
      screened <-
        screened %>%
        filter(!sorter %in% screened_switched_vendor$sorter) %>%
        bind_rows(screened_switched_vendor)
      
      # identify uniprot where commercial Ab is on higher development level than Agrisera
      # these have to later be removed from Agrisera df
      screened_commercial_highest_status_identify <-
        screened %>%
        arrange(desc(assay_order_in_pipeline)) %>%
        group_by(uniprot) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        filter(!str_detect(vendor_a, "agrisera|Agrisera"))
      
      # identify all assays using Agrisera antibodies that were screened and unify their data
      screened_agrisera_prep <-
        screened %>%
        filter(str_detect(vendor_a, "agrisera|Agrisera")) %>% 
        mutate(vendor_a = str_remove(vendor_a, "/abcore")) %>%
        mutate(vendor_a = case_when(vendor_a == "agrisera" ~ paste0("agrisera/", art_no_ag),
                                    .default = vendor_a)) %>%
        mutate(vendor_a = str_remove_all(vendor_a, "agrisera")) %>%
        mutate(vendor_a = str_remove_all(vendor_a, "/")) %>%
        pivot_longer(cols = c(absea, icosagen, bonopus),
                     names_to = "vendor_antigen",
                     values_to = "vendor_category") %>%
        left_join(antigen_categories, by = c("vendor_category" = "antigen_status"))
      
      # from these, select those where the immunogen is different from the immunogen in production by vendors
      # define antigen status but make sure to consider all assays for the same uniprot and their immunogens to determine is there is a new immunogen under development
      screened_agrisera_1 <-
        screened_agrisera_prep %>%
        filter(!is.na(vendor_category)) %>%
        group_by(uniprot) %>%
        mutate(vendor_a_all = paste0(vendor_a, collapse = "; ")) %>%
        filter(vendor_a != vendor_antigen) %>%
        ungroup() %>%
        relocate(vendor_a_all, .after = vendor_a) %>%
        filter(!str_detect(vendor_a_all, vendor_antigen)) %>% 
        arrange(desc(antigen_order_in_pipeline)) %>%
        group_by(uniprot, vendor_a) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        arrange(desc(assay_order_in_pipeline)) %>%
        group_by(uniprot) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        select(names(never_initiated))
      
      # from the rest, select those where immunogen is the same as in production by vendors or where there is no immunogen in production
      screened_agrisera_2 <-
        screened_agrisera_prep %>%
        filter(!uniprot %in% screened_agrisera_1$uniprot) %>%
        arrange(desc(assay_order_in_pipeline)) %>%
        group_by(uniprot) %>%
        slice_head(n = 1) %>%
        mutate(antigen_status_category = "Not initiated") %>%
        select(names(never_initiated))
      
      # combine both agrisera screened parts and remove those where commercial Ab is on higher development level
      screened_agrisera <-
        bind_rows(screened_agrisera_1,
                  screened_agrisera_2) %>%
        filter(!uniprot %in% screened_commercial_highest_status_identify$uniprot)
      
      # select those screened where commercial Ab is on highest level but Agrisera initiative exist in screening, determine antigen status based on result for agrisera assays
      screened_commercial_agrisera_exist <-
        screened_commercial_highest_status_identify %>%
        filter(uniprot %in% screened_agrisera_prep$uniprot) %>%
        left_join(screened_agrisera_1 %>%
                    bind_rows(screened_agrisera_2) %>%
                    select(uniprot, antigen_status_category), by = "uniprot")
      
      # select those screened where commercial Ab is on highest level and no Agrisera initiative exist in screening, determine whether any immunogen is in production
      screened_only_commercial <-
        screened_commercial_highest_status_identify %>%
        filter(!uniprot %in% screened_agrisera$uniprot &
                 !uniprot %in% screened_commercial_agrisera_exist$uniprot) %>%
        pivot_longer(cols = c(absea, icosagen, bonopus),
                     names_to = "vendor_antigen",
                     values_to = "vendor_category") %>%
        left_join(antigen_categories, by = c("vendor_category" = "antigen_status")) %>%
        arrange(desc(antigen_order_in_pipeline)) %>%
        group_by(uniprot) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        mutate(antigen_status_category = if_else(is.na(antigen_status_category), "Never initiated", antigen_status_category)) %>%
        select(names(never_initiated))
      
      # combine all antigen status parts
      status_antigen_final <-
        bind_rows(never_initiated,
                  not_screened_ag_initiated,
                  screened_agrisera,
                  screened_commercial_agrisera_exist,
                  screened_only_commercial)
      
      # combine screening clear status with antigen evaluated status
      status_final <-
        bind_rows(status_screening_clear %>%
                    mutate(antigen_status_category = NA),
                  status_antigen_final) %>%
        mutate(assay_status_category = if_else(assay_status_category == "In product", "Removed from product", assay_status_category)) %>%
        mutate(assay_status_antigen_status = case_when(is.na(antigen_status_category) ~ assay_status_category,
                                                       .default = paste0(assay_status_category, ";", antigen_status_category))) %>%
        mutate(assay_status = if_else(assay_status_category == "In released product", "in product", assay_status))
      
      # bind original data columns for vendor and article nrs from d72
      # create antibody pair type column and identify polyclonal abs
      status_final_extra_columns <-
        status_final %>%
        select(uniprot, assay_status_antigen_status, assay_status_category, assay_status, antigen_status_category, sorter) %>%
        left_join(status_level_adjust %>% select(uniprot, sorter, vendor_a, vendor_b, artnr_a_arm, art_nr_b_arm, art_no_ag, antigen, vendor_ag, issue, comment_status), by = c("uniprot", "sorter")) %>%
        mutate(antibody_pair_type = case_when(is.na(vendor_a) ~ "",
                                              str_detect(vendor_a, "Agrisera| agrisera") & str_detect(vendor_b, "Agrisera| agrisera") ~ "Agrisera",
                                              str_detect(vendor_a, "Agrisera| agrisera")  & !str_detect(vendor_b, "Agrisera| agrisera") ~ "Agrisera + commercial",
                                              !str_detect(vendor_a, "Agrisera| agrisera")  & str_detect(vendor_b, "Agrisera| agrisera") ~ "Agrisera + commercial",
                                              .default = "commercial")) %>%
        mutate(artnr_a_arm_temp = str_remove(artnr_a_arm, "-CUSTOM.*|/CUSTOM.*|_CUSTOM.*"),
               art_nr_b_arm_temp = str_remove(art_nr_b_arm, "-CUSTOM.*|/CUSTOM.*|_CUSTOM.*")) %>%
        mutate(polyclonal = case_when(is.na(vendor_a) ~ "",
                                      artnr_a_arm_temp == art_nr_b_arm_temp & !(str_detect(vendor_a, "Hycult|Hytest|Medix|Abcam")) ~ "yes", # we only get monoclonal antibodies from these vendors, but they can be used on both arms
                                      .default = "")) %>%
        relocate(c(antibody_pair_type, polyclonal), .after = antigen_status_category) %>%
        select(-sorter, -artnr_a_arm_temp, -art_nr_b_arm_temp)
      
      homology <-
        homology_list %>%
        filter(perc_ident >= 80) %>%
        mutate(perc_ident = round(perc_ident, 0)) %>%
        arrange(desc(perc_ident)) %>%
        group_by(query_id, subject_id) %>%
        slice_head(n = 1) %>% # only get the highest homology result for the same uniprot pair
        ungroup() %>%
        arrange(desc(perc_ident)) %>%
        mutate(homolog_percent = paste0(subject_id, ": ", perc_ident, "%")) %>%
        mutate(in_product = if_else(subject_id %in% in_product$uniprot, "homolog_in_product", "homolog_not_in_product")) %>%
        group_by(query_id, in_product) %>%
        summarise(homolog_percent = paste(homolog_percent, collapse = "; ")) %>%
        pivot_wider(names_from = in_product,
                    values_from = homolog_percent) %>%
        select("uniprot" = query_id, homolog_in_product, homolog_not_in_product)
      
      status_final_extra_columns_homology <-
        status_final_extra_columns %>%
        left_join(homology, by = "uniprot") %>%
        relocate(c(homolog_in_product, homolog_not_in_product), .after = antigen_status_category) %>%
        arrange(uniprot)
      
      # export the data
      write_xlsx(status_final_extra_columns_homology, paste0("output/", Sys.Date(), "_assay_status_per_uniprot.xlsx"))
      
    })
    
    cat("✅ The script has run succesfully! Exported results in excel format can be found in the 'output' folder.\n")
    
    return(status_final_extra_columns_homology)
  }
}

