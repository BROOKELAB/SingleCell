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
          "results/test_filter_and_normalize/2020-03-09")
indir <- args[1]          ###indir: input directory. this is a folder output from CellRanger containing the raw matrix (i.e. not the filtered matrix)
out <- args[2]            ###out: output base name (multiple files are output)

##load initial SimpleSingleCell libraries####
library(simpleSingleCell)
library(DropletUtils)
library(scater)
library(EnsDb.Hsapiens.v86)
library(scran)
library(BiocSingular)
library(magrittr)
library(mixtools)
library(ggplot2)


##load object####
sce <- read10xCounts(indir, col.names=TRUE)
dim(sce)
#33702 2211840


##filtering out empty cells####
##call cells: monte carlo p-values; p-value = sig difference from ambient pool rna
##defaultDrops() is alternative that uses 10X method (more conservative, requires cell count)
set.seed(seed)
e.out <- emptyDrops(counts(sce))
#using which() to automatically remove NAs and retain only detected cells
sce <- sce[,which(e.out$FDR <= 0.001)]
dim(sce)
#33702 13598

##calculate preliminary QC stats####
sce <- calculateQCMetrics(sce)  ###no MT control here, cells are infected

##create cell cycle scores####
#should be done before filtering out low-abundance genes
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
cellcycles <- cyclone(sce, pairs=hs.pairs)

sce$CellCycle <- cellcycles$phases
rm(cellcycles)

saveRDS(sce, file = paste0(out, "_sce_postCellCycle.rds"))
sce <- readRDS(paste0(out, "_sce_postCellCycle.rds"))

##get unique rownames####
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)

###correct virus gene NA name:

rownames(sce)[rownames(sce) %in% "23308118"] <- "NA"

##Filter 1: filter cells by min expressed features (i.e. filter matrix columns)####
##remove low feature cell columns
keep.cell <- sce$total_features_by_counts > min.features
sce <- sce[,keep.cell]
remove(keep.cell)
dim(sce)
#33702 13182

##Filter 2: Filter features by min cells (i.e. filter rows)####
##remove low expressed gene rows
keep.feature <- nexprs(sce, byrow=TRUE) >= min.cells
sce <- sce[keep.feature,]
remove(keep.feature)
dim(sce)
#19779 13182

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
#USING ONLY HOST COUNTS

temp1 <- colSums(head(assay(sce, "normcounts"), -8))
sce <- sce[,temp1 < 2*median(temp1)]
rm(temp1)
dim(sce)
#19779 12850

##calculate final QC stats####
sce <- calculateQCMetrics(sce)  ###no MT control here, cells are infected


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
#remove(clusters)
#remove(sf)

##create logged & non-logged normalized values####
sce <- logNormCounts(sce,log = FALSE)
sce <- logNormCounts(sce,log = TRUE)


#### Add in library IDs ####

Lib.ID <- strsplit(sce$Barcode, "-") %>% sapply(function(x) x[2]) %>% factor
levels(Lib.ID) <- c("Bystander","Mock","Infected")
sce$Library <- factor(Lib.ID, levels = c("Mock","Bystander","Infected"))

table(sce$Library, clusters)

#### Calculate total normcount, host normcounts and viral normcounts----

sce$total_normcount <- colSums(assay(sce, "normcounts"))
sce$host_normcount <- colSums(head(assay(sce, "normcounts"), -8))
sce$virus_normcount <- colSums(tail(assay(sce, "normcounts"), 8))

sce$virus_pct <- sce$virus_normcount / sce$total_normcount *100
halfmin_virus <- min(sce$virus_pct[sce$virus_pct>0])
halfmin_virus

sce$virus_pct_log10 <- log10(sce$virus_pct+halfmin_virus)


# 
# 
# x11(height = 10, width = 6)
# layout(matrix(1:3,3,1))
# hist(sce$virus_pct_log10[sce$Library=="Mock"], 1000, xlim = c(-3,2),
#      main = "Mock virus percentage")
# hist(sce$virus_pct_log10[sce$Library=="Bystander"],1000, xlim = c(-3,2),
#      main = "Bystander virus percentage")
# hist(sce$virus_pct_log10[sce$Library=="Infected"],1000, xlim = c(-3,2),
#      main = "Infected virus percentage")

temp <- as.data.frame(colData(sce))

x11()
ggplot(temp, aes(x=total_counts)) + 
  geom_histogram(binwidth=1000) + 
  facet_grid(Library ~ .)
#infected appear to have slightly lower total counts
x11()
ggplot(temp, aes(x=total_normcount)) + 
  geom_histogram(binwidth=1000) + 
  facet_grid(Library ~ .)
#why infected so different after normalization?
#Just due to viral reads?
x11()
ggplot(temp, aes(x=host_normcount)) + 
  geom_histogram(binwidth=1000) + 
  facet_grid(Library ~ .)
#why are regular gene counts shifted so far up?

x11()
ggplot(temp, aes(x=virus_pct)) + 
  geom_histogram() + 
  facet_grid(Library ~ .) 
#most have 0-25% of viral reads

x11()
ggplot(temp, aes(x=virus_pct_log10)) + 
  geom_histogram(binwidth=0.1) + 
  facet_grid(Library ~ .) 

x11()
ggplot(temp, aes(x=virus_pct_log10)) + 
  geom_histogram(binwidth=0.005) + 
  facet_grid(Library ~ .) + 
  ylim(0,40)


#Also check raw counts

sce$host_count <- colSums(head(assay(sce, "counts"), -8))
sce$virus_count <- colSums(tail(assay(sce, "counts"), 8))
sce$virus_pct_raw <- sce$virus_count / sce$total_counts *100


temp <- as.data.frame(colData(sce))

x11()
ggplot(temp, aes(x=host_count)) + 
  geom_histogram(binwidth=1000) + 
  facet_grid(Library ~ .)
#lower, as could be expected

x11()
ggplot(temp, aes(x=virus_pct_raw)) + 
  geom_histogram() + 
  facet_grid(Library ~ .) 
#almost the same?

cor(temp$virus_pct, temp$virus_pct_raw)
#1
all.equal(temp$virus_pct, temp$virus_pct_raw)
#TRUE
#This is probably correct because normalization is within each cell


#Do the libraries drastically differ in the number of non-zero genes?

x11()
ggplot(temp, aes(x=total_features_by_counts)) + 
  geom_histogram() + 
  facet_grid(Library ~ .) 




x11()
ggplot(temp, aes(x=virus_normcount)) + geom_histogram() +
  scale_x_log10()

temp <- data.frame(virus_pct = sce$virus_pct[sce$Library== "Infected"])


x11()
ggplot(temp, aes(x=virus_pct)) + geom_histogram(binwidth=1)

hist(sce$virus_pct[sce$Library=="Infected"], 20)


#look at distributions on log10 scale----
#get half the minimum non-zero value

halfmin_virus <- min(sce$virus_pct[sce$virus_pct>0])
halfmin_virus

sce$virus_pct_log10 <- log10(sce$virus_pct+halfmin_virus)
range(sce$virus_pct_log10)
#-2.821893  1.964439

x11(height = 10, width = 6)
layout(matrix(1:3,3,1))
hist(sce$virus_pct_log10[sce$Library=="Mock"], 1000, xlim = c(-3,2),
     main = "Mock virus percentage")
hist(sce$virus_pct_log10[sce$Library=="Bystander"],1000, xlim = c(-3,2),
     main = "Bystander virus percentage")
hist(sce$virus_pct_log10[sce$Library=="Infected"],1000, xlim = c(-3,2),
     main = "Infected virus percentage")

x11(height = 10, width = 6)
jpeg("results/test_filter_and_normalize/2020-03-12-VirusProportions-AllLibraries.jpeg",
     height = 10, width = 6, units = "in", res = 300, quality = 100)
layout(matrix(1:3,3,1))
hist(sce$virus_pct_log10[sce$Library=="Mock"], 1000, xlim = c(-3,2),
     main = "Mock virus percentage", ylim = c(0,10))
hist(sce$virus_pct_log10[sce$Library=="Bystander"],1000, xlim = c(-3,2),
     main = "Bystander virus percentage", ylim = c(0,10))
hist(sce$virus_pct_log10[sce$Library=="Infected"],1000, xlim = c(-3,2),
     main = "Infected virus percentage")
dev.off()

#what percentage of cells of each type have 0 virus_pct?

table(sce$Library, sce$virus_pct ==0)
#            FALSE TRUE
# Mock        313 4282
# Bystander   234 3375
# Infected   4646    0

table(sce$Library, sce$virus_pct ==0) %>% prop.table(margin = 1)
#                FALSE       TRUE
# Mock      0.06811752 0.93188248
# Bystander 0.06483791 0.93516209
# Infected  1.00000000 0.00000000

min(sce$virus_pct[sce$Library == "Infected"])


#Calculate bimodal index and use + 3 SD of lower peak
#for virus_pct_log10 

mixmdl_virus = normalmixEM(sce$virus_pct_log10[sce$Library=="Infected"],
                           mu = c(-0.25,1), lambda = c(0.2, 0.8))
#NOTE: the above is not determinant and the exact values may
#depend on the number of iterations

x11()
plot(mixmdl_virus,which=2)
lines(density(sce$virus_pct_log10[sce$Library=="Infected"]), lty=2, lwd=2)

mixmdl_virus$lambda
#0.1014708 0.8985292
mixmdl_virus$mu
#-0.3174081  0.9938582
mixmdl_virus$sigma
#0.07165349 0.36666551

i <- which.min(mixmdl_virus$mu)
thresh_infected <- mixmdl_virus$mu[i] + 3*mixmdl_virus$sigma[i]

x11(height = 10, width = 6)
layout(matrix(1:3,3,1))
hist(sce$virus_pct_log10[sce$Library=="Mock"], 1000, xlim = c(-3,2),
     main = "Mock virus percentage", ylim = c(0,10))
abline(v = quantile(sce$virus_pct_log10[sce$virus_pct>0 & sce$Library == "Mock"], 0.95), col = "blue", lty = 2, lwd = 2)
hist(sce$virus_pct_log10[sce$Library=="Bystander"],1000, xlim = c(-3,2),
     main = "Bystander virus percentage", ylim = c(0,10))
abline(v = quantile(sce$virus_pct_log10[sce$virus_pct>0 & sce$Library == "Bystander"], 0.95), col = 2, lty = 2, lwd = 2)
hist(sce$virus_pct_log10[sce$Library=="Infected"],1000, xlim = c(-3,2),
     main = "Infected virus percentage")
abline(v = thresh_infected, col = 3, lty = 2, lwd = 2)


sum(sce$virus_pct_log10[sce$Library=="Bystander"]> 0)

#remove these two obviously infected cells----

sce <- sce[,-which(sce$virus_pct_log10 > 0 & sce$Library=="Bystander")]
dim(sce)
# 19779 12848


#call infected any cell that is > 3 SD----

sce$InfectedStatus <-  ifelse(sce$virus_pct_log10 > thresh_infected, "infected", "notinfected")

table(sce$InfectedStatus, sce$Library)


## Calculate individual virus gene percentages----

n.colData <- ncol(colData(sce)) 
virus_gene_normcounts <- tail(assay(sce, "normcounts"), 8) %>% t() %>% as.matrix()
colData(sce) <- cbind(colData(sce), virus_gene_normcounts/sce$total_normcount * 100)
names(colData(sce))[(n.colData+1):(n.colData+8)] <- paste0(names(colData(sce))[(n.colData+1):(n.colData+8)],"_pct")

#find half the min non-zero value for each----

halfmin_PB2 <- min(sce$PB2_pct[sce$PB2_pct>0])
sce$PB2_pct_log10 <- log10(sce$PB2_pct + halfmin_PB2)
halfmin_PB1 <- min(sce$PB1_pct[sce$PB1_pct>0])
sce$PB1_pct_log10 <- log10(sce$PB1_pct + halfmin_PB1)
halfmin_PA <- min(sce$PA_pct[sce$PA_pct>0])
sce$PA_pct_log10 <- log10(sce$PA_pct + halfmin_PA)
halfmin_HA <- min(sce$HA_pct[sce$HA_pct>0])
sce$HA_pct_log10 <- log10(sce$HA_pct + halfmin_HA)
halfmin_NP <- min(sce$NP_pct[sce$NP_pct>0])
sce$NP_pct_log10 <- log10(sce$NP_pct + halfmin_NP)
halfmin_NA <- min(sce$NA_pct[sce$NA_pct>0])
sce$NA_pct_log10 <- log10(sce$NA_pct + halfmin_NA)
halfmin_M <- min(sce$M_pct[sce$M_pct>0])
sce$M_pct_log10 <- log10(sce$M_pct + halfmin_M)
halfmin_NS <- min(sce$NS_pct[sce$NS_pct>0])
sce$NS_pct_log10 <- log10(sce$NS_pct + halfmin_NS)


x11(width = 6, height = 10)
layout(matrix(1:8, 4, 2, byrow = TRUE))
for (i in grep("log10$", names(colData(sce)))[-1] ){
  hist(colData(sce)[sce$Library == "Infected",i], 1000, xlim = c(-3,2), 
       main = names(colData(sce))[i])
}
#PB2, PB1, and PA have many zeros
#put on same y-axis

x11(width = 6, height = 10)
layout(matrix(1:8, 4, 2, byrow = TRUE))
for (i in grep("log10$", names(colData(sce)))[-1] ){
  hist(colData(sce)[sce$Library == "Infected",i], 1000, xlim = c(-3,2), 
       main = names(colData(sce))[i], ylim = c(0,30))
}



for (i in grep("log10$", names(colData(sce)))[-1] ){
  x11(height = 10, width = 6)
  layout(matrix(1:3,3,1))
  temp1 <- colData(sce)[sce$Library == "Mock",i]
  hist(temp1, 1000, xlim = c(-3,2),
       main = "Mock virus percentage", ylim = c(0,10))
  abline(v = quantile(temp1[temp1>-2.75], 0.95), col = "blue", lty = 2, lwd = 2)
  temp2 <- colData(sce)[sce$Library == "Bystander",i]
  hist(temp2,1000, xlim = c(-3,2),
       main = "Bystander virus percentage", ylim = c(0,10))
  abline(v = quantile(temp2[temp2>-2.75], 0.95), col = 3, lty = 2, lwd = 2)
  hist(colData(sce)[sce$Library == "Infected",i],1000, xlim = c(-3,2),
       main = "Infected virus percentage", ylim = c(0,30))
  abline(v = quantile(temp2[temp2>-2.75], 0.95), col = 3, lty = 2, lwd = 2)
}




  hist(colData(sce)[sce$Library == "Infected",i], 1000, xlim = c(-3,2), 
       main = names(colData(sce))[i], ylim = c(0,30))
}




##Save rds file----

saveRDS(sce, file = paste0(out, "_sce_finalHostVirus.rds"))
#sce <- readRDS(paste0(out, "_sce_finalHostVirus.rds"))

##Create first metadata file ####
meta.1 <- colData(sce)[,c("Barcode","Library","total_features_by_counts","total_counts","CellCycle")]
out.meta <- paste(out,"_1st_metadata.tsv",sep = "")
write.table(meta.1, file = out.meta, sep = "\t")




q()
