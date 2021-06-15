Refactoring load\_data.R
================

-   [Compare “get\_all\_data” function](#compare-get_all_data-function)
-   [Compare get\_expression\_of\_interest and
    get\_expression\_of\_interest\_copy](#compare-get_expression_of_interest-and-get_expression_of_interest_copy)

``` r
knitr::opts_chunk$set()
library(GREAT)
library(ggplot2)
library(data.table)
library(cowplot)
library(ggpubr)
```

    ## 
    ## Attaching package: 'ggpubr'

    ## The following object is masked from 'package:cowplot':
    ## 
    ##     get_legend

``` r
source("../../R/load_data.R")
source("../../R/load_data_copy.R")
```

``` r
# Define data sample path 

id_table_path <- system.file("extdata/sample_data/id_table.RDS", package = "GREAT")
arabidopsis_expression_path <- system.file("extdata/sample_data/arabidopsis_expression.RDS", package = "GREAT")
brassica_rapa_expression_path <- system.file("extdata/sample_data/brassica_rapa_expression.RDS", package = "GREAT")
```

## Compare “get\_all\_data” function

``` r
# data_to_align <- readRDS(arabidopsis_expression_path)
# data_target <- readRDS(brassica_rapa_expression_path)
# id_table <- readRDS(id_table_path)
# colnames_id_table = c("CDS.model", "symbol", "locus_name")
# id_table_unique <- unique(id_table[, ..colnames_id_table]) %>% 
#   dplyr:: mutate_all(.funs = toupper) 

test_get_all_data <- get_all_data(filepath_data_target = brassica_rapa_expression_path, 
                         filepath_data_to_align = arabidopsis_expression_path,
                         filepath_id_table = id_table_path,
                         target_id_table_shared_colname = "CDS.model",
                         target_and_to_align_data_shared_colname = "locus_name",
                         colnames_id_table = c("CDS.model", "symbol", "locus_name"),
                         colnames_wanted = NULL)


test_get_all_data
```

    ##             CDS.model                    sample_id      FPKM accession tissue
    ##   1: BRAA02G018970.3C  ERR_ro18_rna_seq_v3_R18A1_1  0.287365      Ro18   apex
    ##   2: BRAA02G018970.3C  ERR_ro18_rna_seq_v3_R18A1_2  1.069565      Ro18   apex
    ##   3: BRAA02G018970.3C  ERR_ro18_rna_seq_v3_R18A1_3  0.541893      Ro18   apex
    ##   4: BRAA02G018970.3C ERR_ro18_rna_seq_v3_R18A10_1  9.197430      Ro18   apex
    ##   5: BRAA02G018970.3C ERR_ro18_rna_seq_v3_R18A10_2 17.884268      Ro18   apex
    ##  ---                                                                         
    ## 359:        AT5G61850  ERR_ds_klepikova_SRR2073146 29.288324      Col0   apex
    ## 360:        AT5G61850  ERR_ds_klepikova_SRR2073174  4.415715      Col0   apex
    ## 361:        AT5G61850  ERR_ds_klepikova_SRR2073176  2.924128      Col0   apex
    ## 362:        AT5G61850  ERR_ds_klepikova_SRR2073179 12.868319      Col0   apex
    ## 363:        AT5G61850  ERR_ds_klepikova_SRR2106520 12.938941      Col0   apex
    ##      timepoint        dataset est_counts        group   norm.cpm locus_name
    ##   1:        11   ro18_rna_seq         12 Ro18-apex-11  0.3968734  AT1G69120
    ##   2:        11   ro18_rna_seq         47 Ro18-apex-11  1.4147711  AT1G69120
    ##   3:        11   ro18_rna_seq         23 Ro18-apex-11  0.7423984  AT1G69120
    ##   4:        29   ro18_rna_seq        343 Ro18-apex-29 11.3007002  AT1G69120
    ##   5:        29   ro18_rna_seq        749 Ro18-apex-29 23.2055664  AT1G69120
    ##  ---                                                                       
    ## 359:        13 ds_klepikova_4        354 Col0-apex-13 42.9752377  AT5G61850
    ## 360:         9 ds_klepikova_3         40 Col0-apex-09  7.2502672  AT5G61850
    ## 361:         9 ds_klepikova_4         31 Col0-apex-09  5.1475432  AT5G61850
    ## 362:        12 ds_klepikova_3        166 Col0-apex-12 22.2130362  AT5G61850
    ## 363:        12 ds_klepikova_4        177 Col0-apex-12 19.1034466  AT5G61850

``` r
test_get_all_data_copy <- get_all_data_copy(file_path_brassica = brassica_rapa_expression_path, 
                                            file_path_arabidopsis = arabidopsis_expression_path, 
                                            file_path_id_table = id_table_path,
                                            colnames_wanted = NULL)

test_get_all_data_copy
```

    ##             CDS.model                    sample_id      FPKM accession tissue
    ##   1: BRAA02G018970.3C  ERR_ro18_rna_seq_v3_R18A1_1  0.287365      Ro18   apex
    ##   2: BRAA02G018970.3C  ERR_ro18_rna_seq_v3_R18A1_2  1.069565      Ro18   apex
    ##   3: BRAA02G018970.3C  ERR_ro18_rna_seq_v3_R18A1_3  0.541893      Ro18   apex
    ##   4: BRAA02G018970.3C ERR_ro18_rna_seq_v3_R18A10_1  9.197430      Ro18   apex
    ##   5: BRAA02G018970.3C ERR_ro18_rna_seq_v3_R18A10_2 17.884268      Ro18   apex
    ##  ---                                                                         
    ## 359:        AT5G61850  ERR_ds_klepikova_SRR2073146 29.288324      Col0   apex
    ## 360:        AT5G61850  ERR_ds_klepikova_SRR2073174  4.415715      Col0   apex
    ## 361:        AT5G61850  ERR_ds_klepikova_SRR2073176  2.924128      Col0   apex
    ## 362:        AT5G61850  ERR_ds_klepikova_SRR2073179 12.868319      Col0   apex
    ## 363:        AT5G61850  ERR_ds_klepikova_SRR2106520 12.938941      Col0   apex
    ##      timepoint        dataset est_counts        group   norm.cpm locus_name
    ##   1:        11   ro18_rna_seq         12 Ro18-apex-11  0.3968734  AT1G69120
    ##   2:        11   ro18_rna_seq         47 Ro18-apex-11  1.4147711  AT1G69120
    ##   3:        11   ro18_rna_seq         23 Ro18-apex-11  0.7423984  AT1G69120
    ##   4:        29   ro18_rna_seq        343 Ro18-apex-29 11.3007002  AT1G69120
    ##   5:        29   ro18_rna_seq        749 Ro18-apex-29 23.2055664  AT1G69120
    ##  ---                                                                       
    ## 359:        13 ds_klepikova_4        354 Col0-apex-13 42.9752377  AT5G61850
    ## 360:         9 ds_klepikova_3         40 Col0-apex-09  7.2502672  AT5G61850
    ## 361:         9 ds_klepikova_4         31 Col0-apex-09  5.1475432  AT5G61850
    ## 362:        12 ds_klepikova_3        166 Col0-apex-12 22.2130362  AT5G61850
    ## 363:        12 ds_klepikova_4        177 Col0-apex-12 19.1034466  AT5G61850

``` r
setdiff(test_get_all_data, test_get_all_data_copy)
```

    ## Null data.table (0 rows and 0 cols)

## Compare get\_expression\_of\_interest and get\_expression\_of\_interest\_copy

``` r
test_get_expression_of_interest <- get_expression_of_interest(
  filepath_data_target = brassica_rapa_expression_path,
  filepath_data_to_align = arabidopsis_expression_path,
  filepath_id_table = id_table_path,
  target_id_table_shared_colname = "CDS.model",
  target_and_to_align_data_shared_colname = "locus_name",
  colnames_id_table = c("CDS.model", "symbol", "locus_name"),
  colnames_wanted = NULL,
  tissue_wanted = "apex",
  curr_GoIs = c("AT1G69120", "AT5G61850", "AT2G45660"),
  sum_exp_data_target = F, 
  accession_data_to_align = "Col0",
  ids_data_target_colnames = c("CDS.model", "locus_name")
)
```

    ## Note: Using an external vector in selections is ambiguous.
    ## ℹ Use `all_of(target_id_table_shared_colname)` instead of `target_id_table_shared_colname` to silence this message.
    ## ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This message is displayed once per session.

``` r
test_get_expression_of_interest
```

``` r
test_get_expression_of_interest_copy <- get_expression_of_interest_copy(file_path_brassica = brassica_rapa_expression_path, file_path_arabidopsis = arabidopsis_expression_path, file_path_id_table = id_table_path, tissue_wanted = "apex", curr_GoIs =  c("AT1G69120", "AT5G61850", "AT2G45660"), sum_brassicas = F)

test_get_expression_of_interest_copy 
```

``` r
setdiff(test_get_expression_of_interest, test_get_expression_of_interest_copy)
```

    ## Null data.table (0 rows and 0 cols)
