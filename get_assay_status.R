library(janitor)
library(readxl)
library(tidyverse)

uniprot_list <- read_csv("data/2025_03_06_reviewed_uniprot_ids_human.csv")

assay_categories <- read_xlsx("data/status_categories_assay_antigen.xlsx", sheet = "assay") %>%
  clean_names() %>%
  mutate(assay_status = str_to_lower(assay_status))

antigen_categories <- read_xlsx("data/status_categories_assay_antigen.xlsx", sheet = "antigen") %>%
  clean_names() %>%
  mutate(antigen_status = str_to_lower(antigen_status))

in_product <- read_csv("data/2025_10_16_proteins_in_products.csv") %>%
  mutate(in_product = case_when(product == "released_product" ~ "In released product",
                                product == "maximus" ~ "In Maximus",
                                .default = product)) %>%
  select(-product) %>%
  distinct()

# d72 <- read_xlsx("data/DOK084-072 Summary of results in PD084.xlsx", sheet = "All assays", skip = 2) %>%
#   clean_names() %>%
#   rename("assay_status" = status)
# 
# master_list <- read_xlsx("data/Masterlist.xlsx") %>%
#   clean_names() %>%
#   rename("uniprot" = "uniprot_id")

source("functions/get_assay_status_function.R")

master_list <-
  read_xlsx("data/Data driven Priolist with status update_20250911_AEH.xlsx", sheet = "Analysis") %>%
  clean_names()

d72 <- read_xlsx("data/Master_list med 072 status.xlsx", sheet = "Dok_084-72") %>%
  clean_names() %>%
  rename("assay_status" = status)


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
              select(uniprot, "Status AEH"))

write_xlsx(assay_status_resut_compare, "assay_status_resut_compare_2025_10_21.xlsx")

status_compare_count <-
assay_status_resut_compare %>%
  group_by(assay_status_antigen_status, `Status AEH`) %>%
  count()

write_xlsx(status_compare_count, "status_compare_count_2025_10_21.xlsx")

d72_status <- 
  d72 %>%
  select(uniprot, assay_status)

