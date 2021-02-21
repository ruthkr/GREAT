# Creating subsets rds/csv data to test the load functions

# Get original data
rds_rapa <- readRDS("dev/final_data/rds/ro18_chiifu_apex.rds")
id_tab_rapa <- readRDS("dev/reference_data//ID_TABLE_brapa-v3.rds")
klepikova <- readRDS("dev/final_data/rds/klepikova.rds")

# Get information from reference table for specific floral genes
id_table_rapa_subset <- id_tab_rapa %>%
  dplyr::filter(symbol %in% c("SOC1", "AP1", "LFY"))

# Get the brassica accession for those specific floral genes
bra_wanted <- id_table_rapa_subset %>%
  dplyr::pull(CDS.model) %>%
  unique()

# Get the arabidopsis accession for those specific floral genes
ara_wanted <- id_table_rapa_subset %>%
  dplyr::pull(locus_name) %>%
  unique()

# Subsetting the data
klepikova_subset <- klepikova %>%
  dplyr::filter(CDS.model %in% ara_wanted)

rapa_subset <- rds_rapa %>%
  dplyr::filter(CDS.model %in% bra_wanted)


# Save to rds for all
saveRDS(id_table_rapa_subset, "inst/extdata/sample_data/id_table.RDS")
saveRDS(klepikova_subset, "inst/extdata/sample_data/arabidopsis_expression.RDS")
saveRDS(rapa_subset, "inst/extdata/sample_data/brassica_rapa_expression.RDS")
