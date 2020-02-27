#! /usr/bin/env Rscript

##libraries used/installed:
##simpleSingleCell_1.10.1
##EnsDb.Hsapiens.v86_2.99.0  ##edit to different ENSEMBL version if needed
##DropletUtils_1.6.1
##scater_1.14.6
##scran_1.14.5
##BiocSingular_1.2.1

###defaults/settings####
min.features = 400 ##min features/cell for column filter
min.cells = 4      ##min cells/feature for row filter
min.counts = 50    ##min molecule count per cell for column filter OPTIONAL?
perp = 50          ##perplexity for tsne plot
seed = 100

# ##specify which factors to import to object and (optionally) plot in t-sne
# factors2plot <- c("Library", "CellCycle", "TotalVirus", "TotalPB2", "TotalPB1", "TotalPA", "TotalHA", "TotalNP", "TotalNA", "TotalM", "TotalNS", "StatusInfected", "StatusPB2", "StatusPB1", "StatusPA", "StatusHA", "StatusNP", "StatusNA", "StatusM", "StatusNS", "AnyMissingVirusGenes", "AnySingleMissingVirusGene", "NumVirusGenes", "ClusterID", "Doublets")


##get command line####
#args <- commandArgs(TRUE)
## Actually, run interactive for now:
args <- c("results/cbrooke_flu_10xsinglecell3_v2_altRef/flu_singlecell_InfectionStatus_aggr1/outs/raw_gene_bc_matrices_mex/cellranger_hg38_cal0709",
          "results/test_filter_and_normalize/2020-02-26")
indir <- args[1]          ###indir: input directory. this is a folder output from CellRanger containing the raw matrix (i.e. not the filtered matrix)
out <- args[2]            ###out: output base name (multiple files are output)

##load initial SimpleSingleCell libraries####
library(simpleSingleCell)
library(DropletUtils)
library(scater)
library(EnsDb.Hsapiens.v86)
library(scran)
library(BiocSingular)

##load object####
sce <- read10xCounts(indir, col.names=TRUE)

#remove viral genes, which are the last 8 genes####

sce <- sce[1:(nrow(sce)-8),]


##filtering out empty cells####
##call cells: monte carlo p-values; p-value = sig difference from ambient pool rna
##defaultDrops() is alternative that uses 10X method (more conservative, requires cell count)
set.seed(seed)
e.out <- emptyDrops(counts(sce))
#using which() to automatically remove NAs and retain only detected cells
sce <- sce[,which(e.out$FDR <= 0.001)]


##calculate preliminary QC stats####
sce <- calculateQCMetrics(sce)  ###no MT control here, cells are infected

##create cell cycle scores####
#should be done before filtering out low-abundance genes
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
cellcycles <- cyclone(sce, pairs=hs.pairs)

sce$CellCycle <- cellcycles$phases

##get unique rownames####
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)


##Filter 1: filter cells by min expressed features (i.e. filter matrix columns)####
##remove low feature cell columns
keep.cell <- sce$total_features_by_counts > min.features
sce <- sce[,keep.cell]
remove(keep.cell)


##Filter 2: Filter features by min cells (i.e. filter rows)####
##remove low expressed gene rows
keep.feature <- nexprs(sce, byrow=TRUE) >= min.cells
sce <- sce[keep.feature,]
remove(keep.feature)

##Do 1st pass normalization####
##generate quick clusters for normalization
clusters <- quickCluster(sce, use.ranks=FALSE, BSPARAM=IrlbaParam())

##create size factors for normalizing within clusters
sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters)
sf <- sce@int_colData@listData$size_factor
sce <- computeSumFactors(sce, min.mean=0.1,scaling = sf)
remove(clusters)
remove(sf)

##create non-logged normalized values####
sce <- logNormCounts(sce,log = FALSE)

##Call doublets####
#Remove cells with more than twice the median total normalized counts

temp1 <- colSums(assay(sce, "normcounts"))
sce <- sce[,temp1 < 2*median(temp1)]


##add chromosome location info
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID,column="SEQNAME", keytype="GENEID")
rowData(sce)$chr <- location
remove(location)


##Do 2nd pass normalization####
##generate quick clusters for normalization
clusters <- quickCluster(sce, use.ranks=FALSE, BSPARAM=IrlbaParam())

##create size factors for normalizing within clusters
sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters)
sf <- sce@int_colData@listData$size_factor
sce <- computeSumFactors(sce, min.mean=0.1,scaling = sf)
remove(clusters)
remove(sf)

##create non-logged normalized values####
sce <- logNormCounts(sce,log = FALSE)
sce <- logNormCounts(sce,log = TRUE)

##Create first metadata file ####

head(colData(sce))


##Output host-only rds file

out.file <- paste(out,"_HostOnly.rds",sep = "")
saveRDS(sce,file = out.file)



q()
