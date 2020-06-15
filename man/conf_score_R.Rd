\name{conf_score_R}
\alias{conf_score_R}
\title{Calculating the confidence score C for false positive control.}
\usage{
conf_score_R(ref, query, null_expr, gene_num)
}
\arguments{
\item{ref}{The reference dataset. Rows should be cells, columns should be genes.}

\item{query}{The query dataset. Rows should be cells and columns should be genes.The genes are intersection of ref , null and query.}

\item{null_expr}{The control dataset}

\item{gene_num}{The number of common markers of reference set used for false positive control.}
}
\value{
  A vector of confidence scores.
}
\description{
Calculating the confidence score C for false positive control.
}