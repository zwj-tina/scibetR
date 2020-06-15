\name{NullTest_R}
\alias{NullTest_R}
\title{Compute null test.}
\usage{
NullTest_R(ref,query,null_expr,labels,gene_num)
}
\arguments{
\item{ref}{The reference dataset. Rows should be cells, columns should be genes.}

\item{query}{The query dataset. Rows should be cells and columns should be genes.The genes are intersection of ref , null and query.}

\item{null_expr}{The control dataset}

\item{labels}{The reference dataset's labels of each cell.}

\item{gene_num}{The number of common markers of reference set used for false positive control.}
}
\value{
 A vector of null test.
}
\description{
Compute null test.
}