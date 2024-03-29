% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterrnaCrosslink.R
\docType{methods}
\name{clusterrnaCrosslink}
\alias{clusterrnaCrosslink}
\alias{clusterrnaCrosslink,rnaCrosslinkDataSet-method}
\title{clusterrnaCrosslink}
\usage{
clusterrnaCrosslink(cds, cores = 3, stepCount = 2, clusterCutoff = 20)
}
\arguments{
\item{cds}{rnaCrosslinkDataSet object created with rnaCrosslinkDataSet}

\item{cores}{numeric - The number of cores to use}

\item{stepCount}{Stringency for clustering}

\item{clusterCutoff}{The minimum number of reads a cluster requires}
}
\value{
A \code{rnaCrosslinkDataSet} object
}
\description{
This method clusters the duplexes.
}
\examples{
cds = makeExamplernaCrosslinkDataSet()

clusterrnaCrosslink(cds,
                cores = 1,
                stepCount = 1,
                clusterCutoff = 0)
}
