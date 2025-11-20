library(janitor)
library(readxl)
library(tidyverse)

uniprot_list <- read_csv("data/2025_03_06_reviewed_uniprot_ids_human.csv") %>%
  select(uniprot)

assay_categories <- read_xlsx("data/status_categories_assay_antigen.xlsx", sheet = "assay") %>%
  clean_names() %>%
  mutate(assay_status = str_to_lower(assay_status)) %>%
  select(assay_status, assay_synonym_to, assay_order_in_pipeline, assay_status_category, antigen_check_needed)

antigen_categories <- read_xlsx("data/status_categories_assay_antigen.xlsx", sheet = "antigen") %>%
  clean_names() %>%
  mutate(antigen_status = str_to_lower(antigen_status)) %>%
  select(antigen_status, antigen_order_in_pipeline, antigen_status_category)

in_released_product <- read_xlsx("data/AllPublishedAssays-20251024144622.xlsx") %>%
  clean_names() %>%
  rename("uniprot" = "uni_prot") %>%
  select(uniprot)

in_maximus <- read_xlsx("data/Product assay lists_2025-07-21.xlsx", sheet = "Maximus") %>%
  rename("uniprot" = 1) %>%
  select(uniprot) %>%
  distinct()

in_product <-
  in_released_product %>%
  bind_rows(in_maximus) %>%
  distinct() %>%
  mutate(in_product = if_else(uniprot %in% in_released_product$uniprot, "In released product", "In Maximus"))


# d72 <- read_xlsx("data/DOK084-072 Summary of results in PD084.xlsx", sheet = "All assays", skip = 2) %>%
#   clean_names() %>%
#   rename("assay_status" = status) %>%
#   filter(!is.na(uniprot)) %>%
#   filter(!is.na(assay_status)) %>%
#   mutate(across(where(is.character), str_squish))
# 
# master_list <- read_xlsx("data/Masterlist.xlsx") %>%
#   clean_names() %>%
#   rename("uniprot" = "uniprot_id") %>%
#   mutate(across(where(is.character), str_squish))

source("functions/get_assay_status_function.R")

master_list <-
  read_xlsx("data/Data driven Priolist with status update_20250911_AEH.xlsx", sheet = "Analysis") %>%
  clean_names() %>%
  mutate(across(where(is.character), str_squish))

d72 <- read_xlsx("data/Master_list med 072 status.xlsx", sheet = "Dok_084-72") %>%
  clean_names() %>%
  select(uniprot, status, vendor_a, art_no_ag, artnr_a_arm, art_nr_b_arm, vendor_b, issue, antigen, vendor_ag, comment_status, sorter) %>%
  rename("assay_status" = status) %>%
  filter(!is.na(uniprot)) %>%
  filter(!is.na(assay_status)) %>%
  mutate(across(where(is.character), str_squish))

assay_status_resut <- get_assay_status(
  uniprot_list,
  assay_categories,
  antigen_categories,
  in_product,
  d72,
  master_list
)

assay_status_resut %>%
  select(assay_status_antigen_status) %>%
  distinct() %>%
  View()

assay_status_resut %>%
  group_by(uniprot) %>%
  count() %>%
  arrange(desc(n))

angelica <- read_csv("data/Data driven Priolist with status update_To Kerhan_20250912(DA Priolist ).csv")

assay_status_resut_compare <-
  assay_status_resut %>%
  left_join(angelica %>% 
              select(uniprot, "Status AEH")) %>%
  relocate(`Status AEH`, .after = uniprot)

write_xlsx(assay_status_resut_compare, "assay_status_resut_compare_2025_11_17.xlsx")

status_compare_count <-
  assay_status_resut_compare %>%
  group_by(assay_status_antigen_status, `Status AEH`) %>%
  count() %>%
  arrange(desc(n))

write_xlsx(status_compare_count, "status_compare_count_2025_11_17.xlsx")

