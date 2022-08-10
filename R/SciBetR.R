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
#' @name Gambler_R
#' @param test_r Rows should be cells and columns should be genes.
#' @param prob_r  trained scibet model  Rows should be genes and columns should be cell types.The genes must be matched with test_r
#' @param ret_tab   return list or matrix
#' @return  'cellType'  or matrix
Gambler_R <- function(test_r, prob_r,ret_tab=FALSE){
  total <- as.matrix(test_r) %*% as.matrix(prob_r)
  cellType <- c()
  for(i in 1:nrow(total)){
    index <- which.max(total[i,])
    cellType <- c(cellType,colnames(total)[index])
  }
  mm <- 2 ^ (total-max(total))
  out <- mm/ sum(mm)
  if(ret_tab){
    return(out)
  }
  return(cellType)
}
#' Train SciBet model and generate a "Bet" function for cell type prediction.
#' @name Learn_R
#' @usage Learn_R(expr, geneset, k=1000, a=5)
#' @param expr The reference expression dataframeThe expression dataframe, with rows being cells, and columns being genes. The last column should be "label".
#' @param geneset A user-defined set of genes for prediction.
#' @param k See [SelectGene()] for details.
#' @return A function that takes two inputs: `test` and `result`. See [Bet()] for details.
#' @export

Learn_R <- function(expr,geneset=NULL,k=1000){
  labels <- factor(expr$label)
  if(is.null(geneset)){
    geneset <- SelectGene_R(expr,k)
  }
  labell <- levels(labels)
  expr_select <- expr[,geneset]
  #average
  label_total <- matrix(0,length(geneset),1)
  for(i in labell){
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
  genes <- rownames(prob)
  function(test, result="list"){
    have_genes <- intersect(genes,colnames(test))
    testa <- log1p(as.matrix(test[,have_genes])) / log(2)
    switch(result,
           list = Gambler_R(testa, prob[have_genes, ], FALSE),
           table = {
             out <- Gambler_R(testa, prob[have_genes, ],TRUE)
             rownames(out) <- have_genes
        return(out)
        }
      )
    }
  }
#' Classify cells of a given query dataset using a reference dataset.
#' @description SciBet main function. Train SciBet with the reference dataset to assign cell types for the query dataset.
#' @name SciBet_R
#' @usage SciBet_R(train, test, k=1000, result=c("list", "table"))
#' @param train The reference dataset, with rows being cells, and columns being genes. The last column should be "label".
#' @param test The query dataset. Rows should be cells and columns should be genes.
#' @param k Number of genes to select according to the total entropy differneces among cell types.
#' @examples SciBet_R(train.matr, query.matr)
#' @export
SciBet_R <- function(train, test, k=1000, result="list"){
  Learn_R(train, NULL, k)(test, result)
}

#' Compute expression entropy.
#' @name Entropy_R
#' @usage Entropy_R(expr, window=120, low = 2000)
#' @param expr The expression dataframe. Rows should be cells and columns should be genes.
#' @param window The window size for expression value discretization.
#' @param low The lower limit for normalizing expression entropy
#' @return A dataframe..
#' @export

Entropy_R <- function(expr, window=120, low = 2000){
  expr <- as.data.frame(expr)

  ent_res <- tibble(
    gene = colnames(expr),
    mean.expr = colMeans(expr)
  ) %>%
    dplyr::filter(mean.expr < 6000)

  expr <- expr[,ent_res$gene]
  #
  #out <- GenEntr(expr, window, n_threads)
  #
  n_cell <- nrow(expr)
  n_gene <- ncol(expr)
  out <- c()
  for(i in 1:n_gene){
    states <- ceiling(1000000.0/window)
    discretize <- vector(mode="numeric", length=states)
    li <- ceiling(expr[,i]/window)+1
    number <- table(li)
    for(j in names(number)){
      a <- number[j]
      discretize[as.numeric(j)] <- a
    }
    sumAll <- sum(discretize)
    discretizeN <- discretize[discretize!=0]
    out <- c(out,-sum(discretizeN/sumAll * log(discretizeN/sumAll)))
  }

  ent_res %>%
    dplyr::mutate(entropy = out) %>%
    dplyr::mutate(fit = 0.18*log(0.03*mean.expr + 1)) -> ent_res

  ent_res %>%
    dplyr::filter(mean.expr > low) %>%
    dplyr::mutate(k = entropy/fit) %>%    #linear normalization of expression entropy
    dplyr::pull(k) %>%
    quantile(0.75) %>%
    as.numeric() -> k

  ent_res <- ent_res %>% dplyr::mutate(norm_ent = entropy/k)

  return(ent_res)
}

#' Calculating the nulltest.
#' @name NullTest_R
#' @usage NullTest_R(ref, query, null_expr, labels,gene_num = 500)
#' @param ref The reference dataset. Rows should be cells, columns should be genes.
#' @param query The query dataset. Rows should be cells and columns should be genes.The genes are intersection of ref , null and query.
#' @param labels The reference dataset's labels of each cell.
#' @param gene_num The number of common markers of reference set used for false positive control.
#' @return A vector of null test.
#' @export
NullTest_R <- function(ref,query,null_expr,labels,gene_num){
   n_ref <- nrow(ref)
   n_query <- nrow(query)
   n_gene <- ncol(ref)
   labell <- levels(labels)
   n_label <- length(labell)
   label_total <- matrix(0,ncol(ref),1)
   for(i in labell){
     label_TPM <- ref[labels==i,]
     label_mean <- colMeans(label_TPM)
     a <- matrix(label_mean,,1)
     colnames(a) <- i
     label_total <- cbind(label_total,a)
   }
   label_total <- label_total[,-1]
   rownames(label_total)<- colnames(ref)
   e1 <- c()
   ds <- c()
   e1 =  rowSums(label_total)/n_label
   ds = log(e1+1)-log(null_expr+1)
   sum_e1 <- sum(e1)
   sum_null <- sum(null_expr)
   e1 <- log((e1+1)/(sum_e1+n_gene))-log((null_expr+1)/(sum_null+n_gene))
   dss <- sort(ds,decreasing = TRUE)
   prob <- vector(mode = "numeric",length = n_query)
   for(j in 1:n_query){
       prob[j] <- sum(e1[names(dss)[1:gene_num]]*query[j,names(dss)[1:gene_num]])
   }
   return(prob)
}
#' Calculating the confidence score C for false positive control.
#' @name conf_score_R
#' @usage conf_score(ref, query, null_expr, gene_num = 500)
#' @param ref The reference dataset. Rows should be cells, columns should be genes and last column should be "label".
#' @param query The query dataset. Rows should be cells and columns should be genes.
#' @param gene_num The number of common markers of reference set used for false positive control.
#' @return A vector of confidence scores.
#' @export
conf_score_R <- function(ref, query, null_expr, gene_num){
  labels <- factor(ref$label)
  genes <- Reduce(intersect, list(colnames(ref), colnames(query), names(null_expr)))
  ref <- log1p(as.matrix(ref[, genes])) / log(2)
  query <- log1p(as.matrix(query[, genes])) / log(2)
  a <- NullTest_R(ref, query, null_expr[genes], labels, gene_num)
  b <- NullTest_R(ref, ref, null_expr[genes], labels, gene_num)
  prob <- a/max(b)
  prob[prob > 1] <- 1
  return(prob)
}

#' Make predictions with the trained SciBet model.
#' @name Bet_R
#' @usage Bet(expr, result=c("list", "table"))
#' @param expr The expression matrix or dataframe for prediction. Rows should be cells and columns should be genes.
#' @param result Return a "list" of predicted labels, or a "table" of probabilities of each tested cell belonging to each label.
#' @return A vector or a dataframe for the classification result.
#' @export
Bet_R <- function(expr, result="list"){
}

#' Export model as matrix.
#' @details If you don't plan to export the models for usage on other platforms, simply saving the Bet function in a RData file would also work in R.
#' @name ExportModel_R
#' @usage ExportModel(Bet)
#' @param Bet A prediction function generated by SciBet.
#' @return A matrix.
#' @export
ExportModel_R <- function(Bet_R){
  prob <- environment(Bet_R)$prob
  rownames(prob) <- environment(Bet_R)$genes
  colnames(prob) <- environment(Bet_R)$labell
  return(prob)
}
#' Generate Bet function from a model matrix.
#' @name LoadModel_R
#' @usage LoadModel(x, genes, labels)
#' @param x A SciBet model in the format of a matrix.
#' @param genes (Optional).
#' @param labels (Optional).
#' @return A Bet function.
#' @export
LoadModel_R <- function(x, genes=NULL, labels=NULL){
  prob <- x
  if (is.null(genes))
    genes <- rownames(x)
  if (is.null(labels))
    labels <- colnames(x)
  function(expr, result="list"){
    have_genes <- intersect(genes,colnames(expr))
    expra <- log1p(as.matrix(expr[,have_genes])) / log(2)
    switch(result,
           list = Gambler_R(expra, prob[have_genes,],FALSE),
           table = {
             out <- Gambler_R(expra, prob[have_genes, ], TRUE)
             rownames(out) <- have_genes
             return(out)
           }
    )
  }
}


#' Process scibet.core
#' @name pro.core
#' @usage pro.core(scibet.core)
#' @param scibet.core A SciBet core
#' @return A processed SciBet core
#' @export
pro.core <- function(scibet.core){
  cell.type <- unname(unlist(scibet.core[,1]))
  scibet.core <- as.data.frame(t(scibet.core[,-1]))
  colnames(scibet.core) <- cell.type
  return(as.matrix(scibet.core))
}


#' Expression patterns of informative genes across cell types.
#' @name Marker_heatmap
#' @usage Marker_heatmap(expr, gene)
#' @param expr The expression dataframe. Rows should be cells, columns should be genes and last column should be "label".
#' @param gene A vector of informative genes.
#' @return A figure.
#' @export
Marker_heatmap <- function(expr, gene){
  expr <- expr[,c(gene,'label')]
  type_expr <- expr %>%
    tidyr::nest(-label) %>%
    dplyr::rename(expr = data) %>%
    dplyr::mutate(colmeans = purrr::map(
      .x = expr,
      .f = function(.x){colMeans(.x)}))

  type_expr$colmeans %>%
    as.data.frame() %>%
    tibble::remove_rownames() %>%
    t() %>%
    as.data.frame() %>%
    tibble::remove_rownames() -> type_mean_expr

  rownames(type_mean_expr) <- type_expr$label
  colnames(type_mean_expr) <- colnames(expr)[-ncol(expr)]

  sub_expr <- type_mean_expr
  sub_expr <- sub_expr
    as_tibble() %>%
    dplyr::mutate_all(funs((. - mean(.))/sd(.))) %>%
    t()
  colnames(sub_expr) <- type_expr$label
  get_label <- function(num){
    v <- sub_expr[num,]
    colnames(sub_expr)[which(v == max(v))]
  }
  sub_expr <- sub_expr %>%
    tibble::as.tibble() %>%
    dplyr::mutate(group = purrr::map_chr(1:length(gene), get_label))
  sub_expr <- as.data.frame(sub_expr)
  rownames(sub_expr) <- gene
  sub_expr <- sub_expr %>%
    dplyr::mutate(gene = gene) %>%
    tidyr::gather(key = 'cell_type', value = 'zscore', -group, -gene) %>%
    dplyr::arrange(group, desc(zscore))
  sub_expr %>%
    ggplot(aes(factor(gene, levels = unique(gene)),
               factor(cell_type, levels = sort(unique(cell_type), decreasing = T)))) +
    geom_point(aes(size = zscore, colour = zscore)) +
    theme(
      strip.text.x = element_blank(),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 13),
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(color="black", angle = -90, hjust = 0),
      panel.background = element_rect(colour = "black", fill = "white"),
      panel.grid = element_line(colour = "grey", linetype = "dashed"),
      panel.grid.major = element_line(
        colour = "grey",
        linetype = "dashed",
        size = 0.2
      )
    ) +
    facet_grid(. ~ group, scales = "free", space = "free") +
    scale_colour_distiller(palette = "RdYlBu") +
    labs(
      x = '',
      y = ''
    ) -> p

  return(p)
}

#' Heatmap of classification result.
#' @name Confusion_heatmap
#' @usage Confusion_heatmap(ori, prd)
#' @param ori A vector of the original labels for each cell in the test set.
#' @param prd A vector of the predicted labels for each cell in the test set.
#' @return A heatmap for the confusion matrix of the classification result.
#' @export
Confusion_heatmap <- function(ori, prd){
  tibble(
    ori = ori,
    prd = prd
  ) %>%
    dplyr::count(ori, prd) %>%
    tidyr::spread(key = prd, value = n) -> cross.validation.filt

  cross.validation.filt[is.na(cross.validation.filt)] = 0
  cross.validation.filt[,-1] <- round(cross.validation.filt[,-1]/rowSums(cross.validation.filt[,-1]),2)
  cross.validation.filt <- cross.validation.filt %>%
    tidyr::gather(key = 'prd', value = 'Prob', -ori)

  cross.validation.filt %>%
    ggplot(aes(ori,prd,fill = Prob)) +
    geom_tile() +
    theme(axis.title = element_text(size = 0)) +
    theme(axis.text = element_text(size = 10)) +
    theme(legend.title = element_text(size = 0)) +
    theme(legend.text = element_text(size = 10)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    theme(axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black",angle = 45, hjust = 1)) +
    scale_fill_viridis() -> p

  return(p)
}


#' Heatmap for the confusion matrix of the classification with the false postive control.
#' @name Confusion_heatmap_negctrl
#' @usage Confusion_heatmap_negctrl(res, cutoff = 0.4)
#' @param res Classification result.
#' @param cutoff The cutoff of confifence score C.
#' @return A heatmap for the confusion matrix of the classification result with the false postive control.
#' @export
Confusion_heatmap_negctrl <- function(res, cutoff = 0.4){
  res %>%
    dplyr::mutate(prd = ifelse(c_score < cutoff, 'unassigned', prd)) %>%
    dplyr::count(ori, prd) %>%
    tidyr::spread(key = prd, value = n) -> cla.res

  cla.res[is.na(cla.res)] = 0
  cla.res[,-1] <- round(cla.res[,-1]/rowSums(cla.res[,-1]),2)
  cla.res <- cla.res %>% tidyr::gather(key = 'prd', value = 'Prob', -ori)
  label <- cla.res$ori %>% unique()
  cla.res %>%
    ggplot(aes(prd, factor(ori, levels = c(label[-3],'Neg.cell')), fill = Prob)) +
    geom_tile(colour = 'white', lwd = 0.5) +
    theme(axis.title = element_text(size = 12)) +
    theme(axis.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 0)) +
    theme(legend.text = element_text(size = 12)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    theme(axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black", angle = 50, hjust = 1)) +
    scale_fill_material('blue')
}


