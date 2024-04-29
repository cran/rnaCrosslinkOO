## ----include = FALSE----------------------------------------------------------
library(rmarkdown)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)




## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo = T, results = 'hide',message=FALSE---------------------------------
#install.packages("rnaCrosslinkOO")
# Load the rnaCrosslink-OO Library 
library(rnaCrosslinkOO)

## ----echo = T, results = 'hide',message=FALSE---------------------------------
# Here are the other libraries on which rnaCrosslinkOO relies
#library(seqinr)
#library(GenomicRanges)
#library(ggplot2)
#library(reshape2)
#library(MASS)
#library(ggplot2)
#library(doParallel)
#library(igraph)
#library(R4RNA)
#library(RColorBrewer)
#library(heatmap3)
#library(mixtools)
#library(TopDom)
library(tidyverse)
#library(RRNA)
#library(ggrepel)

## -----------------------------------------------------------------------------
# Set up the sample table
sampleTableRow1 = c(system.file("extdata", 
                                "s1.txt", 
                                package="rnaCrosslinkOO"), "s", "1", "s1")
sampleTableRow2 = c(system.file("extdata", 
                                "c1.txt", 
                                package="rnaCrosslinkOO"), "c", "1", "c1")
sampleTable2 = rbind.data.frame(sampleTableRow1, sampleTableRow2)

# add the column names 
colnames(sampleTable2) = c("file", "group", "sample", "sampleName")

sampleTable2

## -----------------------------------------------------------------------------
rna = c("ENSG000000XXXXX_NR003286-2_RN18S1_rRNA")

## -----------------------------------------------------------------------------
path18SFata <- system.file("extdata", 
                           "18S.fasta", 
                           package="rnaCrosslinkOO")

rnaRefs = list()
rnaRefs[[rna]] = read.fasta(path18SFata)


## -----------------------------------------------------------------------------
path18SFata <- system.file("extdata", 
                           "ribovision18S.txt", 
                           package="rnaCrosslinkOO")
known18S = read.table(path18SFata,
                      header = F)

## -----------------------------------------------------------------------------
pathShape <- system.file("extdata",
                         "reactivities.txt", 
                         package="rnaCrosslinkOO")
shape = read.table(pathShape,
                      header = F)

## ----eval=FALSE---------------------------------------------------------------
#  rnaCrosslinkQC(sampleTable2,
#                 directory = ".",
#                 topTranscripts = F)
#  

## -----------------------------------------------------------------------------
# load the object
cds = rnaCrosslinkDataSet(rnas = rna,
                      sampleTable = sampleTable2,
                      subset = "all",
                      sample = "all")
# be aware there are extra options here
# including:
# subset - allows you to choose specific read sizes (this affects resolution and accuracy)
# sample - Choose the same ammount of reads for each sample (useful for comparisons)

## -----------------------------------------------------------------------------
# Check status of instance 
cds

## -----------------------------------------------------------------------------
# Returns the size of the RNA
rnaSize(cds)

## -----------------------------------------------------------------------------
# Returns the sample table 
sampleTable(cds)

## -----------------------------------------------------------------------------
# Returns indexes of the samples in the control and not control groups
group(cds)

## -----------------------------------------------------------------------------
# Get the sample names of the instance
sampleNames(cds)

## ----eval=FALSE,echo=TRUE-----------------------------------------------------
#  # Return the InputFiles slot
#  InputFiles(cds)

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  # Return the matrixList slot
#  matrixList(cds)

## -----------------------------------------------------------------------------
data = getData(x    = cds,              # The object      
               data =  "InputFiles",       # The Type of data to return     
               type = "original")[[1]] # The stage of the analysis for the return data

head(data)

## -----------------------------------------------------------------------------
# Returns the RNAs with highest number of assigned reads 
# regardless of whether it is an Inter or Intra - RNA interaction. 
topTranscripts(cds,
               2)  # The number of entries to return


## -----------------------------------------------------------------------------
# Returns the RNAs that interact with the RNA of interest
topInteracters(cds, # The rnaCrosslinkDataSet instance
               1)   # The number of entries to return

## -----------------------------------------------------------------------------
# Returns the Interacions with the highest number of assigned reads
topInteractions(cds, # The rnaCrosslinkDataSet instance
                2)   # The number of entries to return


## -----------------------------------------------------------------------------
features = featureInfo(cds) # The rnaCrosslinkDataSet instance


# Counts for features at the transcript level
features$transcript

# Counts for features at the family level (last field with  "_" delimited IDs)
features$family


## -----------------------------------------------------------------------------
# Cluster the reads
clusteredCds = clusterrnaCrosslink(cds = cds,     # The rnaCrosslinkDataSet instance 
                               cores = 1,         # The number of cores
                               stepCount = 2,     # The number of steps in the random walk
                               clusterCutoff = 2) # The minimum number of reads for a cluster to be considered

## -----------------------------------------------------------------------------
# Check status of instance 
clusteredCds

## -----------------------------------------------------------------------------
# Returns the number of clusters in each sample
clusterNumbers(clusteredCds)

## -----------------------------------------------------------------------------
# Returns the number reads in clusters
readNumbers( clusteredCds)

## -----------------------------------------------------------------------------
getData(clusteredCds,        # The object             
        "clusterTableList",  # The Type of data to return     
        "original")[[1]]     # The stage of the analysis for the return data

## ----eval=FALSE---------------------------------------------------------------
#  getData(clusteredCds,         # The object
#          "clusterGrangesList", # The Type of data to return
#          "original")[[1]]      # The stage of the analysis for the return data

## -----------------------------------------------------------------------------
# Trim the Clusters
trimmedClusters = trimClusters(clusteredCds = clusteredCds, # The rnaCrosslinkDataSet instance 
                               trimFactor = 1,              # The cutoff for cluster trimming (see above)
                               clusterCutoff = 0)          # The minimum number of reads for a cluster to be considered

## -----------------------------------------------------------------------------
# Check status of instance 
trimmedClusters

## -----------------------------------------------------------------------------
# Returns the number of clusters in each sample
clusterNumbers(trimmedClusters)

## -----------------------------------------------------------------------------
# Returns the number reads in clusters
readNumbers( trimmedClusters)

## -----------------------------------------------------------------------------
#plotClusterAgreement(trimmedClusters,
#                     "trimmedClusters")

## -----------------------------------------------------------------------------
#plotClusterAgreementHeat(trimmedClusters,
#                         "trimmedClusters")

## ----error=FALSE,eval = FALSE, results=FALSE,message=FALSE--------------------
#  # Fold the RNA in part of whole
#  foldedCds = foldrnaCrosslink(trimmedClusters,
#                           rnaRefs = rnaRefs,
#                           start = 1600,
#                           end = 1869,
#                           shape = 0,
#                           ensembl = 20,
#                           constraintNumber  = 30,
#                           evCutoff = 5)

## ----eval = FALSE-------------------------------------------------------------
#  # Check status of instance
#  foldedCds

## -----------------------------------------------------------------------------
# Plot heatmaps for each sample
#plotMatrices(cds = cds,         # The rnaCrosslinkDataSet instance 
#             type = "original", # The "analysis stage"
#             directory = 0,     # The directory for output (0 for standard out)
#             a = 1,             # Start coord for x-axis
#             b = rnaSize(cds),  # End coord for x-axis
#             c = 1,             # Start coord for y-axis
#             d = rnaSize(cds),  # End coord for y-axis
#             h = 5)             # The height of the image (if saved)


## -----------------------------------------------------------------------------
plotCombinedMatrix(cds,
           type1 = "original",
           type2 = "noHost",
           b = rnaSize(cds),
           d = rnaSize(cds))

## -----------------------------------------------------------------------------
trimmedClusters

## -----------------------------------------------------------------------------
# Plot heatmaps for all samples combined and all controls combined
plotMatricesAverage(cds = trimmedClusters, # The rnaCrosslinkDataSet instance 
             type1 = "trimmedClusters", # The "analysis stage" to plot on the upper half of the heatmap
             type2 = "original",        # The "analysis stage" to plot on the lower half of the heatmap
             directory = 0,     # The directory for output (0 for standard out)
             a = 1,             # Start coord for x-axis
             b = rnaSize(cds),  # End coord for x-axis
             c = 1,             # Start coord for y-axis
             d = rnaSize(cds),  # End coord for y-axis
             h = 5)             # The height of the image (if saved)

## ----eval = FALSE-------------------------------------------------------------
#  domainDF = data.frame()
#  for(j in c(20,30,40,50,60,70)){
#      #for(i in which(sampleTable(cds)$group == "s")){
#  
#      timeMats = as.matrix(getData(x = cds,
#                                   data = "matrixList",
#                                   type = "noHost")[[1]])
#  
#      timeMats = timeMats/ (sum(timeMats)/1000000)
#      tmp = tempfile()
#      write.table(timeMats, file = tmp,quote = F,row.names = F, col.names = F)
#  
#      tdData2 = readHiC(
#          file = tmp,
#          chr = "rna18s",
#          binSize = 10,
#          debug = getOption("TopDom.debug", FALSE)
#      )
#  
#      tdData =  TopDom(
#          tdData2 ,
#          window.size = j,
#          outFile = NULL,
#          statFilter = TRUE,
#          debug = getOption("TopDom.debug", FALSE)
#      )
#  
#      td = tdData$domain
#      td$sample = sampleTable(cds)$sampleName[1]
#      td$window = j
#      domainDF = rbind.data.frame(td, domainDF)
#  
#  }
#  
#  
#  
#  ggplot(domainDF) +
#      geom_segment(aes(x = from.coord/10,
#                       xend = to.coord/10, y = as.factor(sub("s","",sample)),
#                       yend = (as.factor(sub("s","",sample)) ), colour = tag),
#                   size  = 20, alpha = 0.8) +
#      facet_grid(window~.)+
#      theme_bw()
#  

## ----eval = FALSE-------------------------------------------------------------
#  plotEnsemblePCA(foldedCds,
#                  labels = T, # plot labels for structures
#                  split = T)  # split samples over different facets (T/f)
#  

## ----eval = FALSE-------------------------------------------------------------
#  plotComparisonArc(foldedCds = foldedCds,
#                    s1 = "c1",            # The sample of the 1st structure
#                    s2 = "s1",            # The sample of the 2nd structure
#                    n1 = 13,               # The number of the 1st structure
#                    n2 = 16)               # The number of the 2nd structure

## ----eval = F-----------------------------------------------------------------
#  plotStructure(foldedCds = foldedCds,
#                rnaRefs = rnaRefs,
#                s = "s1",          # The sample of the structure
#                n = 1)             # The number of the structure

## -----------------------------------------------------------------------------

getInteractions(cds,
                "ENSG00000XXXXXX_NR003287-2_RN28S1_rRNA") %>%
    mutate(sample =sub("\\d$","",sample) )%>%
    group_by(rna,Position,sample)%>%
    summarise(sum =  sum(depth)) %>%
    ggplot()+
    geom_area(aes(x = Position,
                  y = sum, 
                  fill = sample), 
              stat = "identity")+
    facet_grid(sample~.) +
    theme_bw()


## -----------------------------------------------------------------------------
getReverseInteractions(cds,
                       rna) %>%
    mutate(sample =sub("\\d$","",sample) )%>%
    group_by(rna,Position,sample)%>%
    summarise(sum =  sum(depth)) %>%
    ggplot()+
    geom_area(aes(x = Position,
                  y = sum, 
                  fill = sample), 
                    stat = "identity")+
    facet_grid(sample~.)+
    theme_bw()

## -----------------------------------------------------------------------------
plotInteractions(cds,
                 rna = "ENSG000000XXXXX_NR003286-2_RN18S1_rRNA",
                 interactor = "ENSG00000XXXXXX_NR003287-2_RN28S1_rRNA",
                 b = "max",
                 d = "max")

## -----------------------------------------------------------------------------
plotInteractionsAverage(cds,
                 rna = "ENSG000000XXXXX_NR003286-2_RN18S1_rRNA",
                 interactor = "ENSG00000XXXXXX_NR003287-2_RN28S1_rRNA",
                 b = "max",
                 d = "max")

## -----------------------------------------------------------------------------
expansionSize = 5
knownMat = matrix(0, nrow = rnaSize(cds), ncol = rnaSize(cds))
for(i in 1:nrow(known18S)){
    knownMat[ (known18S$V1[i]-expansionSize):(known18S$V1[i]+expansionSize),
              (known18S$V2[i]-expansionSize):(known18S$V2[i]+expansionSize)] =
        knownMat[(known18S$V1[i]-expansionSize):(known18S$V1[i]+expansionSize),
                 (known18S$V2[i]-expansionSize):(known18S$V2[i]+expansionSize)] +1
}
knownMat = knownMat + t(knownMat)


## -----------------------------------------------------------------------------
# use compare known to gett he known and not know clusters
knowClusteredCds = compareKnown(trimmedClusters, # The rnaCrosslinkDataSet instance 
                                knownMat, # A contact matrix of know interactions
                                "trimmedClusters") # The analysis stage of clustering to compare 

knowClusteredCds

## -----------------------------------------------------------------------------
# Plot heatmaps for all samples combined and all controls combined
plotMatricesAverage(cds = knowClusteredCds, # The rnaCrosslinkDataSet instance 
             type1 = "KnownAndNovel", # The "analysis stage"
             directory = 0,     # The directory for output (0 for standard out)
             a = 1,             # Start coord for x-axis
             b = rnaSize(cds),  # End coord for x-axis
             c = 1,             # Start coord for y-axis
             d = rnaSize(cds),  # End coord for y-axis
             h = 5)             # The hight of the image (if saved)


## -----------------------------------------------------------------------------
# Get the number of clusters for each analysis Stage
clusterNumbers(knowClusteredCds)

## -----------------------------------------------------------------------------
# Get the number of reads in each cluster for each analysis stage
readNumbers(knowClusteredCds)

## ----eval = FALSE-------------------------------------------------------------
#  head(compareKnownStructures(foldedCds,
#                              known18S)) # the comarison set

## ----eval = FALSE-------------------------------------------------------------
#  ggplot(compareKnownStructures(foldedCds, known18S)) +
#      geom_hline(yintercept = c(0.5,0.25,0.75,0,1),
#                 colour = "grey",
#                 alpha = 0.2)+
#      geom_vline(xintercept = c(0.5,0.25,0.75,0,1),
#                 colour = "grey",
#                 alpha = 0.2)+
#      geom_point(aes(x = sensitivity,
#                     y = precision,
#                     size = as.numeric(as.character(unlist(foldedCds@dgs))),
#                     colour = str_sub(structureID,
#                                      start = 1 ,
#                                      end = 2))) +
#      xlim(0,1)+
#      ylim(0,1)+
#      theme_classic()
#  

## ----eval = FALSE-------------------------------------------------------------
#  # make an example dataset
#  cds = makeExamplernaCrosslinkDataSet()
#  cds

## ----eval = FALSE-------------------------------------------------------------
#  # Have a look at the reads for each sample
#  plotMatricesAverage(cds = cds,
#                      type1 = "original",
#                      directory = 0,
#                      a = 1,
#                      b = rnaSize(cds),
#                      c = 1,
#                      d = rnaSize(cds),
#                      h = 5)
#  

## ----eval = FALSE-------------------------------------------------------------
#  topInteracters(cds,2)

## ----eval = FALSE-------------------------------------------------------------
#  featureInfo(cds)

## ----eval = FALSE-------------------------------------------------------------
#  clusteredCds = clusterrnaCrosslink(cds = cds,
#                                 cores = 1,
#                                 stepCount = 2,
#                                 clusterCutoff = 0)
#  
#  
#  
#  
#  
#  
#  trimmedClusters = trimClusters(clusteredCds = clusteredCds,
#                                 trimFactor = 1,
#                                 clusterCutoff = 0)

## ----eval = FALSE-------------------------------------------------------------
#  plotMatricesAverage(cds = clusteredCds,
#                      type1 = "originalClusters",
#                      directory = 0,
#                      a = 1,
#                      b = rnaSize(cds),
#                      c = 1,
#                      d = rnaSize(cds),
#                      h = 5)

## ----eval = FALSE-------------------------------------------------------------
#  plotMatricesAverage(cds = trimmedClusters,
#                      type1 = "trimmedClusters",
#                      directory = 0,
#                      a = 1,
#                      b = rnaSize(cds),
#                      c = 1,
#                      d = rnaSize(cds),
#                      h = 5)
#  

## ----eval = FALSE-------------------------------------------------------------
#  
#  fasta = paste(c(rep('A',25),
#                  rep('T',25),
#                  rep('A',10),
#                  rep('T',23)),collapse = "")
#  
#  header = '>transcript1'
#  
#  
#  fastaFile = tempfile()
#  writeLines(paste(header,fasta,sep = "\n"),con = fastaFile)
#  
#  
#  rnaRefs = list()
#  rnaRefs[[rnas(cds)]] = read.fasta(fastaFile)
#  rnaRefs
#  
#  
#  
#  foldedCds = foldrnaCrosslink(trimmedClusters,
#                           rnaRefs = rnaRefs,
#                           start = 1,
#                           end = 83,
#                           shape = 0,
#                           ensembl = 5,
#                           constraintNumber  = 1,
#                           evCutoff = 1)
#  

## ----eval = FALSE-------------------------------------------------------------
#  # make an example table of "know" interactions
#  file = data.frame(V1 = c(6),
#                    V2 = c(80))
#  
#  compareKnownStructures(foldedCds,file)
#  

## ----eval = FALSE-------------------------------------------------------------
#  
#  plotEnsemblePCA(foldedCds, labels = T,split = T)

## ----eval = FALSE-------------------------------------------------------------
#  plotComparisonArc(foldedCds = foldedCds,"s1","c1",2,3)

## ----eval = FALSE-------------------------------------------------------------
#  plotStructure(foldedCds = foldedCds,"c1",1)
#  
#  

