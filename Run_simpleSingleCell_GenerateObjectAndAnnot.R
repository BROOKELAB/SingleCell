#! /usr/bin/env Rscript

##libraries used/installed:
##simpleSingleCell_1.8.0
##EnsDb.Hsapiens.v86_2.99.0  ##edit to different ENSEMBL version if needed
##DropletUtils_1.4.1
##scater_1.12.2
##scran_1.12.1
##BiocSingular_1.0.0

###defaults/settings
min.features = 400 ##min features/cell for column filter
min.cells = 4      ##min cells/feature for row filter
min.counts = 50    ##min molecule count per cell for column filter OPTIONAL?
perp = 50          ##perplexity for tsne plot
do.plot1 = 0       ##generate first plot (top 10 bio-variance genes (non-virus))
do.plot2 = 0       ##generate second plot(s) (t-sne for various factors)
seed = 100

##specify which factors to import to object and (optionally) plot in t-sne
factors2plot <- c("Library", "CellCycle", "TotalVirus", "TotalPB2", "TotalPB1", "TotalPA", "TotalHA", "TotalNP", "TotalNA", "TotalM", "TotalNS", "StatusInfected", "StatusPB2", "StatusPB1", "StatusPA", "StatusHA", "StatusNP", "StatusNA", "StatusM", "StatusNS", "AnyMissingVirusGenes", "AnySingleMissingVirusGene", "NumVirusGenes", "ClusterID", "Doublets")


##get command line
args <- commandArgs(TRUE)
indir <- args[1]          ###indir: input directory. this is a folder output from CellRanger containing the raw matrix (i.e. not the filtered matrix)
out <- args[2]            ###out: output base name (multiple files are output)
metafile <- args[3]       ###metafile: this is currently a tab-delimited table with cell ID rows and metadata columns to be added to the SCE object
factorfile <- args[4]     ###factorfile: file containing single column of factor headers from metatable to import to object

##load initial SimpleSingleCell libraries
library(simpleSingleCell)
library(DropletUtils)
library(scater)
library(EnsDb.Hsapiens.v86)
library(scran)
library(BiocSingular)
library(magrittr)


#sce <- readRDS("results/test_GenerateObject/2019-09-26-sce.rds")


##load object
sce <- read10xCounts(indir, col.names=TRUE)


##add chromosome location info
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, column="SEQNAME", keytype="GENEID")
rowData(sce)$chr <- location
remove(location)


##load metadata if present
if (file.exists(metafile)){
  lib <- read.table(metafile,header=TRUE,sep="\t")
  if (file.exists(factorfile)){
    factors2plot <- scan(factorfile,what = "character")
  }
  ##add metadata to singlecellexperiment object
  colData(sce) <- DataFrame(lib)
  remove(lib)
} else {
  cellIDs <- colnames(sce)
  out.file <- paste(out,"_AllCellIDs.tsv",sep = "")
  write.table(cellIDs,file = out.file,sep = "\t")
  remove(cellIDs)
}

##filtering out empty cells
##call cells: monte carlo p-values; p-value = sig difference from ambient pool rna
##defaultDrops() is alternative that uses 10X method (more conservative, requires cell count)
set.seed(seed)
e.out <- emptyDrops(counts(sce))
#using which() to automatically remove NAs and retain only detected cells
sce <- sce[,which(e.out$FDR <= 0.001)]

##calculate preliminary QC stats
sce <- calculateQCMetrics(sce)  ###no MT control here, cells are infected


##filter 1: filter out doublets if column is present
if (file.exists(metafile)){
  keep.doublets <- sce$Doublets == 0
  sce <- sce[,which(keep.doublets)]
  remove(keep.doublets)
}


##filter 2: filter cells by min expressed features (i.e. filter matrix columns)
##remove low feature cell columns
keep.cell <- sce$total_features_by_counts > min.features
sce <- sce[,keep.cell]
remove(keep.cell)


## get sample info

sce$Sample <- substring(sce$Barcode, 18) %>% as.numeric()
sce$Library <- c("Bystander","Mock","Infected")[sce$Sample]



##create cell cycle scores

hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
cellcycles <- cyclone(sce, pairs=hs.pairs)
sce$phase <- cellcycles$phases
rm(cellcycles)


##get unique rownames
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)

## fix NA gene name

rownames(sce)[nrow(sce)-2] <- "NA"





##Filter Option 2: Filter features by min cells (i.e. filter rows)
##remove low expressed gene rows
keep.feature <- nexprs(sce, byrow=TRUE) >= min.cells
sce <- sce[keep.feature,]
remove(keep.feature)

##generate quick clusters for normalization
clusters <- quickCluster(sce, use.ranks=FALSE, BSPARAM=IrlbaParam())

##create size factors for normalizing within clusters
sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters)
sf <- sce@int_colData@listData$size_factor
sce <- computeSumFactors(sce, min.mean=0.1,scaling = sf)
remove(clusters)
remove(sf)

##create normalized log-expression values
sce <- normalize(sce)
sce <- normalize(sce,return_log = FALSE)


## get porportion values from normcounts

normcounts <- as.matrix(assay(sce,"normcounts"))

assay(sce,"props") <- t(t(normcounts) /  colSums(normcounts)) 

rm(normcounts)


## Calculate proportion of all reads from virus and call infection status

virusPropTotal <- colSums(tail(as.matrix(assay(sce,"props")),8))

#use 90th pct of bystander proportion as threshold for "infected"

By.no0.90pctile <- quantile(virusPropTotal[virusPropTotal > 0 & sce$Library == "Bystander"], probs = 0.9)

sce$InfectedStatus <-  ifelse(virusPropTotal > By.no0.90pctile, "infected", "notinfected")


## Calculate individual virus gene proportions and call present/absent

virusPropGenes <- tail(as.matrix(assay(sce,"props")),8) %>% t()

Vgenes.info <- matrix(NA, nrow = 8, ncol = 2, dimnames = list(colnames(virusPropGenes), c("NumNon0By","pctile90")))

for (i in 1:8){
  temp <- virusPropGenes[virusPropGenes[ ,i] > 0 & sce$Library == "Bystander" ,i]
  Vgenes.info[i,1] <- length(temp)
  Vgenes.info[i,2] <- quantile(temp, probs = 0.9)
}

Vgenes.status <- data.frame(matrix(NA, nrow = ncol(sce), ncol = 8, dimnames = list(colnames(sce), paste0(colnames(virusPropGenes), "_Status"))))

for (i in 1:8) {
  Vgenes.status[,i] <- ifelse(virusPropGenes[,i] > Vgenes.info[i,2], "Present", "Absent")
}


## Add in AnyMissingGenes	AnySingleMissingGene	NumVirusGenes

Vgenes.status$AnyMissingGenes <- dplyr::case_when(rowSums(Vgenes.status[,1:8] =="Absent") == 8 ~ "notInfected",
                                                  rowSums(Vgenes.status[,1:8] =="Absent") == 0 ~ "notMissing",
                                                  rowSums(Vgenes.status[,1:8] =="Absent") > 0 & rowSums(Vgenes.status[,1:8] =="Absent") < 0 ~ "Missing")
                                                  
Vgenes.status$AnySingleMissingGene <- dplyr::case_when(rowSums(Vgenes.status[,1:8] =="Absent") == 8 ~ "notInfected",
                                                  rowSums(Vgenes.status[,1:8] =="Absent") == 7 ~ "One",
                                                  rowSums(Vgenes.status[,1:8] =="Absent") < 7 ~ "notOne")


Vgenes.status$NumVirusGenes <- rowSums(Vgenes.status[,1:8] =="Present")



options(stringsAsFactors = FALSE)
colData(sce) <- cbind(colData(sce), Vgenes.status)
options(stringsAsFactors = TRUE)


## calculate doublets based on total host norm counts (minus virus)

hostNormTotal <- colSums(head(as.matrix(assay(sce,"normcounts")),-8))

#separately by library or all together?!?



if (file.exists(metafile)){
  ##save Seurat-safe object
  out.file <- paste(out,".rds",sep = "")
  saveRDS(sce,file = out.file)
} else {
  normcounts <- as.matrix(assay(sce,"normcounts"))
  out.file <- paste(out,"_FilteredMatrix.tsv",sep = "")
  write.table(normcounts,file = out.file, sep = "\t")
  cellIDs <- colnames(sce)
  out.file <- paste(out,"_FilteredCellIDs.tsv",sep = "")
  write.table(cellIDs,file = out.file,sep = "\t")
  remove(cellIDs)
}


if (do.plot1){
  ##model technical noise, assumes noise is Poisson
  new.trend <- makeTechTrend(x=sce)
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


##make Scran t-sne plots
if (do.plot2){
  ##choosing dimension number in data
  sce <- denoisePCA(sce, technical=new.trend, BSPARAM=IrlbaParam())
  
  ##t-sne
  sce <- runTSNE(sce, use_dimred="PCA", perplexity=perp, rand_seed=seed)
  
  ##clustering on PCAs
  snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
  clusters <- igraph::cluster_walktrap(snn.gr)
  sce$Cluster <- factor(clusters$membership)
  for (i in 1:length(factors2plot)){
    plot2 <- plotTSNE(sce, colour_by=factors2plot[i])
    out.plot2 <- paste(out,"_",factors2plot[i],".pdf",sept = "")
    pdf(out.plot2)
    plot(plot2)
    dev.off()
    remove(plot2)
  }
}

q()
