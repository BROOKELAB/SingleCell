#! /usr/bin/env Rscript

##get command line
args <- commandArgs(TRUE)
indir <- args[1]          ###indir: input directory. this is a folder output from CellRanger containing the raw matrix (i.e. not the filtered matrix)
metafile <- args[2]       ###metafile: this is currently "VirusCategories_forR_NewCategories_wCC.txt"
out <- args[3]            ###out: output base name (multiple files are output)

###defaults/settings
min.features = 400 ##min features/cell for column filter
min.cells = 4      ##min cells/feature for row filter
min.counts = 50    ##min 
perp = 50          ##perplexity for tsne plot
do.plot1 = 0       ##generate first plot (top 10 bio-variance genes (non-virus))
do.plot2 = 1       ##generate second plot(s) (t-sne for various factors)
##specify which factors to plot in t-sne
factors2plot <- c("library", "cellcycle", "TotalVirus", "TotalPB2", "TotalPB1", "TotalPA", "TotalHA", "TotalNP", "TotalNA", "TotalM2", "TotalM1", "TotalNS1", "TotalNEP", "StatusInfected", "StatusPB2", "StatusPB1", "StatusPA", "StatusHA", "StatusNP", "StatusNA", "StatusM2", "StatusM1", "StatusNS1", "StatusNEP", "MissingVirusGenes", "SingleMissingVirusGene", "NumMissingVirusGenes", "Cluster")

##load first libraries
library(simpleSingleCell)

##load object
sce <- read10xCounts(indir, col.names=TRUE)

##load cellcycle and other metadata
lib <- read.table(metafile,header=TRUE,sep="\t")

##add metadata to singlecellexperiment object
colData(sce)$library <- lib$library
colData(sce)$cellcycle <- lib$cellcycle
colData(sce)$TotalVirus <- lib$Virus_proportion
colData(sce)$TotalPB2 <- lib$PB2_proportion
colData(sce)$TotalPB1 <- lib$PB1_proportion
colData(sce)$TotalPA <- lib$PA_proportion
colData(sce)$TotalHA <- lib$HA_proportion
colData(sce)$TotalNP <- lib$NP_proportion
colData(sce)$TotalNA <- lib$NA_proportion
colData(sce)$TotalM2 <- lib$M2_proportion
colData(sce)$TotalM1 <- lib$M1_proportion
colData(sce)$TotalNS1 <- lib$NS1_proportion
colData(sce)$TotalNEP <- lib$NEP_proportion
colData(sce)$StatusInfected <- lib$infected_cell
colData(sce)$StatusPB2 <- lib$PB2_status
colData(sce)$StatusPB1 <- lib$PB1_status
colData(sce)$StatusPA <- lib$PA_status
colData(sce)$StatusHA <- lib$HA_status
colData(sce)$StatusNP <- lib$NP_status
colData(sce)$StatusNA <- lib$NA_status
colData(sce)$StatusM2 <- lib$M2_status
colData(sce)$StatusM1 <- lib$M1_status
colData(sce)$StatusNS1 <- lib$NS1_status
colData(sce)$StatusNEP <- lib$NEP_status
colData(sce)$MissingVirusGenes <- lib$Infected_MissingGenes
colData(sce)$SingleMissingVirusGene <- lib$Infected_SingleMissingGene
colData(sce)$NumMissingVirusGenes <- lib$Infected_NumMissingGenes
remove(lib)

##get unique rownames
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)

##add chromosome location info
library(EnsDb.Hsapiens.v86)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID,column="SEQNAME", keytype="GENEID")
rowData(sce)$chr <- location
remove(location)

##filtering out empty cells
##call cells: monte carlo p-values; p-value = sig difference from ambient pool rna
##defaultDrops() is alternative that uses 10X method (more conservative, requires cell count)
set.seed(100)
e.out <- emptyDrops(counts(sce))
#using which() to automatically remove NAs and retain only detected cells
sce <- sce[,which(e.out$FDR <= 0.01)]

##calculate preliminary QC stats
sce <- calculateQCMetrics(sce)  ###no MT control here

##filter option 1: filter cells by min expressed features (i.e. filter columns)
##remove low feature cell columns
keep.cell <- sce$total_features_by_counts > min.features
sce <- sce[,keep.cell]
remove(keep.cell)

##Filter Option 2: Filter features by min cells (i.e. filter rows)
##remove low expressed gene rows
keep.feature <- nexprs(sce, byrow=TRUE) >= min.cells
sce <- sce[keep.feature,]
remove(keep.feature)

##generate quick clusters for normalization
clusters <- quickCluster(sce, method="igraph", min.mean=0.1,irlba.args=list(maxit=1000))
##create size factors for normalizing within clusters
sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters)
remove(clusters)
##create normalized log-expression values
sce <- normalize(sce)

##model technical noise, assumes noise is Poisson
new.trend <- makeTechTrend(x=sce)
if (do.plot1){
  ##generate top (non-virus) 10 gene biological variance plots
  fit <- trendVar(sce, use.spikes=FALSE, loess.args=list(span=0.05))
  fit0 <- fit
  fit$trend <- new.trend
  dec <- decomposeVar(fit=fit)
  top.dec <- dec[order(dec$bio, decreasing=TRUE),]
  plot1 <- plotExpression(sce, features=rownames(top.dec)[11:20])
  ##print out plot
  out.plot1 <- paste(out,"_TopGeneBioVariance.pdf",sept = "")
  pdf(out.plot1)
  plot(plot1)
  dev.off()
  remove(plot1)
}


##choosing dimension number in data
sce <- denoisePCA(sce, technical=new.trend, approx=TRUE)

##t-sne
sce <- runTSNE(sce, use_dimred="PCA", perplexity=perp, rand_seed=100)

##clustering on PCAs
snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
clusters <- igraph::cluster_walktrap(snn.gr)
sce$Cluster <- factor(clusters$membership)

##make t-sne plots
if (do.plot2){
  for (i in 1:length(factors2plot)){
    plot2 <- plotTSNE(sce, colour_by=factors2plot[i])
    out.plot2 <- paste(out,"_",factors2plot[i],".pdf",sept = "")
    pdf(out.plot2)
    plot(plot2)
    dev.off()
    remove(plot2)
  }
}

##save simpleSingleCell object
out.file <- paste(out,".rds",sep = "")
saveRDS(sce,file = out.file)

q()
