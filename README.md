# scibetR
A portable and fast single cell type identifier

## Installation Guide
**Installing dependency package**  
Before installing SciBet, the dependency packages should be installed first:
```
install.packages("Rcpp")
install.packages("RcppEigen")
install.packages("ggsci")
install.packages("viridis")
install.packages("tidyverse")
```
**Installing SciBet**  
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
1. For reference set, rows should be cells, column should be genes and the last column should be "label" (TPM). 2. For query set, rows should be cells and column should be genes (TPM).
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

prd <- SciBet(train_set, test_set[,-ncol(test_set)])
```
