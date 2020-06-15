\name{Entropy_R}
\alias{Entropy_R}
\title{Compute expression entropy.}
\usage{
Entropy_R(expr, window=120, low = 2000)
}
\arguments{
\item{expr}{The expression dataframe. Rows should be cells and columns should be genes.}

\item{window}{The window size for expression value discretization.}

\item{low}{The lower limit for normalizing expression entropy}

\item{result}{Return "table" or "list".}
}
\value{
 A dataframe including mean.expr,entropy,fit and norm_ent
}
\description{
Compute expression entropy.
}
