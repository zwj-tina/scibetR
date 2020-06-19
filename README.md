# scibetR
Pure R version of **scibet**, a portable and fast single cell type identifier

## Installation Guide

**Installing scibetR**  
To install scibetR, run:
```
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("zwj-tina/scibetR")
```

## Tutorial
### Library
```
library(ggplot2)
library(tidyverse)
library(scibetR)
library(viridis)
library(ggsci)
```
### Load the data
For expression matrix (TPM), rows should be cells and the last column should be "label".
```
expr <- 
```
### E(ntropy)-test for supervised gene selection
```
etest_gene <- SelectGene_R(expr, k = 50)
etest_gene
etest_gene
```
### scibetR: Single Cell Identifier Based on Entropy Test
1. For reference set, rows should be cells, column should be genes and the last column should be "label" (TPM).
2. For query set, rows should be cells and column should be genes (TPM).
example:
```
tibble(
  ID = 1:nrow(expr),
  label = expr$label
) %>%
  dplyr::sample_frac(0.7) %>%
  dplyr::pull(ID) -> ID

train_set <- expr[ID,]      #construct reference set
test_set <- expr[-ID,]      #construct query set

prd <- SciBet_R(train_set, test_set[,-ncol(test_set)])
```

### False positive control
Due to the incomplete nature of reference scRNA-seq data collection, cell types excluded from the reference dataset may be falsely predicted to be a known cell type. By applying a null dataset as background, SciBet controls the potential false positives while maintaining high prediction accuracy for cells with types covered by the reference dataset (positive cells).
For the purposes of this example, these three datasets are used to get started.
```
null <- readr::read_rds('~/null.rds.gz')
reference <- readr::read_rds('~/reference.rds.gz')
query <- readr::read_rds('~/query.rds.gz')
```

For query set, “negative cells” account for more than 60%.
```
ori_label <- query$label
table(ori_label)
```

The confidence score of each query cell is calculated with the function conf_score_R.
```
query <- query[,-ncol(query)]
c_score <- conf_score_R(ref = reference, query = query, null_expr = null, gene_num = 500)
```


