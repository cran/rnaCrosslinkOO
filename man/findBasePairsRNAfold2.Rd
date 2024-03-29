% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/foldrnaCrosslinkMethodsAndHelpers.R
\docType{methods}
\name{findBasePairsRNAfold2}
\alias{findBasePairsRNAfold2}
\title{findBasePairsRNAfold2}
\usage{
findBasePairsRNAfold2(startPos, endPos, seqs, fasta)
}
\arguments{
\item{startPos}{Start of the cluster side x}

\item{endPos}{End of the cluster side x}

\item{seqs}{Sequence of x}

\item{fasta}{\code{rnaRefs}}
}
\value{
A table of clusters and coordinates with folds
}
\description{
Folds the clusters using Vienna RNA duplex
}
