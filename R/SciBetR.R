#' @param expr Rows should be cells and the last column should be "label".
#' @param num_top Number of genes to select according to the total entropy differneces among cell types.
#' @return  'out'  selected gene list(class "character)
SelectGene_R <- function(expr, k = 1000) {
  labels <- factor(expr$label)
  labels_set <- levels(labels)
  #average 等同于aggregate函数，但aggregate耗时太长
  label_total <- matrix(0,ncol(expr)-1,1)
  for(i in labels_set){
    label_TPM <- expr[expr$label==i,][,-ncol(expr)]
    label_mean <- colMeans(label_TPM)
    a <- matrix(label_mean,,1)
    colnames(a) <- i
    label_total <- cbind(label_total,a)
  }
  label_total <- label_total[,-1]
  #E-test
  log_E <- log2(rowMeans(label_total+1))
  E_log <- rowMeans(log2(label_total+1))
  t_scores <- log_E - E_log
  out <- cbind(label_total,t_scores)
  rownames(out) <-colnames(expr)[-ncol(expr)]
  out <- out[order(-out[,"t_scores"]),]
  select <- out[1:k,]
  out <- rownames(select)
  return(out)
}
#' @param expr Rows should be cells and the last column should be "label".
#' @param geneset  selected gene
#'  @param k  selected gene number
#' @return  'prob'  scibet matrix Rows are genes and the columns are genes(class matrix)
Learn_R <- function(expr,geneset=NULL,k=1000){
  labels <- factor(expr$label)
  if(is.null(geneset)){
    geneset <- SelectGene_R(expr,k)
  }
  labels_set <- levels(labels)
  expr_select <- expr[,geneset]
  #average
  label_total <- matrix(0,length(geneset),1)
  for(i in labels_set){
    label_TPM <- expr_select[labels==i,]
    label_mean <- colSums(log2(label_TPM+1))
    a <- matrix(label_mean,,1)
    colnames(a) <- i
    label_total <- cbind(label_total,a)
  }
  label_total <- label_total[,-1]
  rownames(label_total) <- geneset
  #calculate probability and log transform
  label_t <- t(label_total)
  prob <- log2(label_t+1)-log2(rowSums(label_t)+length(geneset))
  prob <- t(prob)
  return(prob)
}
#' @param train Rows should be cells and the last column should be "label".
#' @param test Rows should be cells and columns should be genes.
#'  @return  'cellType' the predicted cell types (class character)
SciBet_R <- function(train,test,k=1000){
  scibet_train <- Learn_R(train,NULL,k)
  genes_train <- rownames(scibet_train)
  genes_test <- colnames(test)
  gene_index <- intersect(genes_train,genes_test)
  test_matrix <- log2(test[,gene_index] + 1)
  train_matrix <- scibet_train[gene_index,]
  total <- as.matrix(test_matrix) %*% as.matrix(train_matrix)
  cellType <- c()
  for(i in 1:nrow(total)){
    index <- which.max(total[i,])
    cellType <- c(cellType,colnames(total)[index])
  }
  return(cellType)
}

