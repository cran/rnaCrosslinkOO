% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterrnaCrosslinkMethodsAndHelpers.R
\docType{methods}
\name{compareKnown}
\alias{compareKnown}
\alias{compareKnown,rnaCrosslinkDataSet-method}
\title{compareKnown}
\usage{
compareKnown(trimmedClusters, knownMat, type)
}
\arguments{
\item{trimmedClusters}{a \code{rnaCrosslinkDataSet} object, 
run \code{trimClusters()} on the  \code{rnaCrosslinkDataSet} first}

\item{knownMat}{Matrix - A marix(ncol = lengthRNA,nrow = lengthRNA) where a
value in matrix[x,y] would indicate a known interation between nucleotide 
x and nucleotide y}

\item{type}{string - the Analysis stage of clusters you would like to compare you can find 
available types by just running the objects name}
}
\value{
Returns a \code{rnaCrosslinkClusteredDataSet} object

The 3 attributes matrixList, clusterTableList and clusterGrangesList 
will gain the \code{types} "known" and "novel" and "knownAndNovel"
}
\description{
This method compares the current object to a know structure.run 
\code{trimClusters()} on the  \code{rnaCrosslinkDataSet} first
}
\examples{
cds = makeExamplernaCrosslinkDataSet()

clusteredCds = clusterrnaCrosslink(cds,
                cores = 1,
                stepCount = 1,
                clusterCutoff = 0)
knownMat = matrix(0, ncol = rnaSize(cds), nrow = rnaSize(cds))
knownMat[7,27] = 1
# use compare known to gett he known and not know clusters
knowClusteredCds = compareKnown(clusteredCds,
                                knownMat,
                                "original")
clusterNumbers(knowClusteredCds)                 
                

}
