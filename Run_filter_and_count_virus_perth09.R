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
args <- c("results/perth_1609/cbrooke_flu_10xsinglecell_4_altRef/flu_singlecell_InfectionStatus_aggr4/outs/raw_feature_bc_matrix/",
          "results/test_filter_and_normalize/2020-03-10-perth09")
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


###Fix NA gene ####
#The gene name NA is preventing the data from being read in
#Manually change it:

dir(indir)
temp <- read.delim(paste0(indir, "features.tsv"), header = FALSE, as.is = TRUE)
tail(temp)
temp[33700,1] <- "NA_vir"
temp[33700,2] <- "NA_vir"

write.table(temp, file = paste0(indir, "features.tsv"), quote = FALSE,
            row.names = FALSE, col.names = FALSE, sep = "\t")

#manually gzip



##load object####
sce <- read10xCounts(indir, col.names=TRUE)
dim(sce)
#33702 20384640


#keep viral genes, which are the last 8 genes####
#sce <- sce[1:(nrow(sce)-8),]


##filtering out empty cells####
##call cells: monte carlo p-values; p-value = sig difference from ambient pool rna
##defaultDrops() is alternative that uses 10X method (more conservative, requires cell count)
set.seed(seed)
e.out <- emptyDrops(counts(sce))
#using which() to automatically remove NAs and retain only detected cells
sce <- sce[,which(e.out$FDR <= 0.001)]
dim(sce)
#33702 10803

##calculate preliminary QC stats####
sce <- calculateQCMetrics(sce)  ###no MT control here, cells are infected

##create cell cycle scores####
#should be done before filtering out low-abundance genes
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
cellcycles <- cyclone(sce, pairs=hs.pairs)

sce$CellCycle <- cellcycles$phases
rm(cellcycles)

saveRDS(sce, file = paste0(out, "_sce_postCellCycle.rds"))


#RESTART here 1 ----
sce <- readRDS(paste0(out, "_sce_postCellCycle.rds"))

#### Add in library IDs ####

Lib.ID <- strsplit(sce$Barcode, "-") %>% sapply(function(x) x[2]) %>% factor
levels(Lib.ID) <- c("Bystander","Mock","Infected")
sce$Library <- factor(Lib.ID, levels = c("Mock","Bystander","Infected"))
table(sce$Library)
#Mock Bystander  Infected 
#4862      3502      2439

##get unique rownames####
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)

###correct virus gene NA name:

rownames(sce)[rownames(sce) %in% "NA_vir"] <- "NA"

grep("^IFN", rownames(sce), value = TRUE)
#  [1] "IFNLR1"   "IFNGR1"   "IFNB1"    "IFNW1"    "IFNA21"   "IFNA4"   
#  [7] "IFNA7"    "IFNA10"   "IFNA16"   "IFNA17"   "IFNA14"   "IFNA5"   
# [13] "IFNA6"    "IFNA13"   "IFNA2"    "IFNA8"    "IFNA1"    "IFNE"    
# [19] "IFNK"     "IFNG-AS1" "IFNG"     "IFNL3"    "IFNL2"    "IFNL1"   
# [25] "IFNAR2"   "IFNAR1"   "IFNGR2" 



##Filter 1: filter cells by min expressed features (i.e. filter matrix columns)####
##remove low feature cell columns
keep.cell <- sce$total_features_by_counts > min.features
sce <- sce[,keep.cell]
remove(keep.cell)
dim(sce)
#33702 10059
table(sce$Library)
#Mock Bystander  Infected 
#4495      3166      2398


##Filter 2: Filter features by min cells (i.e. filter rows)####
##remove low expressed gene rows
keep.feature <- nexprs(sce, byrow=TRUE) >= min.cells
sce <- sce[keep.feature,]
remove(keep.feature)
dim(sce)
#19235 10059

grep("^IFN", rownames(sce), value = TRUE)
#"IFNLR1" "IFNGR1" "IFNL2"  "IFNL1"  "IFNAR2" "IFNAR1" "IFNGR2"


# Create graph for doublet calling

sce$doublet <- colSums(head(assay(sce, "counts"), -8))

library(ggridges)
library(Seurat)
so <- CreateSeuratObject(counts = as.matrix(assay(sce, "counts")), 
                         project = "cal07", meta.data = as.data.frame(colData(sce))[,-1])

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

ggplotColours(n = 3)
# "#F8766D" "#00BA38" "#619CFF"

so$Library <- factor(so$Library, levels = c("Infected","Bystander","Mock"))

table(so$Library)
# Infected Bystander      Mock 
#     2398      3166      4495 

x11(5,5)
RidgePlot(so, features = "doublet", group.by = "Library", assay = "RNA",
          same.y.lims = TRUE, cols = c("#00BA38","#F8766D","#619CFF"),
          log = FALSE, slot = "data") +  
  geom_vline(xintercept = 2*median(so$doublet), size = 1) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  xlab("UMI counts of host transcripts") + ylab("") + NoLegend() + ggtitle("Perth09") +
  theme(plot.title = element_text(hjust = 1))

ggsave("results/test_filter_and_normalize/2020-05-04-perth09_doublet_calling.pdf")

rm(so)


temp1 <- colSums(head(assay(sce, "counts"), -8))
x11(height = 10, width = 6)
layout(matrix(1:3,3,1))
hist(temp1[sce$Library=="Mock"], 1000, xlim = c(750,96000),
     main = "Mock")
abline(v = 2*median(temp1), col =2)
hist(temp1[sce$Library=="Bystander"],1000, xlim = c(750,96000),
     main = "Bystander")
abline(v = 2*median(temp1), col =2)
hist(temp1[sce$Library=="Infected"],1000, xlim = c(750,96000),
     main = "Infected")
abline(v = 2*median(temp1), col =2)


##Filter 3: Call doublets####
#Remove cells with more than twice the median total counts
#USING ONLY HOST COUNTS

#temp1 <- colSums(head(assay(sce, "normcounts"), -8))
temp1 <- colSums(head(assay(sce, "counts"), -8))
sce <- sce[,temp1 < 2*median(temp1)]
rm(temp1)
dim(sce)
#19235  9329
table(sce$Library)
#Mock Bystander  Infected 
#4140      2802      2387


##calculate final QC stats####
sce <- calculateQCMetrics(sce)  ###no MT control here, cells are infected


##add chromosome location info
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID,column="SEQNAME", keytype="GENEID")
rowData(sce)$chr <- location
remove(location)


#### Save copy for import into Seurat ####

out <- "results/test_filter_and_normalize/2020-03-30-perth09"
saveRDS(sce, file = paste0(out, "_sce_noNorm.rds"))


#RESTART here 2 ----
sce <- readRDS(paste0(out, "_sce_noNorm.rds"))


##DO NOT DO normalization####

#### Calculate host counts and viral counts ----

sce$host_counts <- colSums(head(assay(sce, "counts"), -8))
sce$virus_counts <- colSums(tail(assay(sce, "counts"), 8))
sce$virus_pct <- sce$virus_counts / sce$total_counts *100


table(sce$Library, sce$virus_pct >0) %>% prop.table(margin = 1)
#                  FALSE         TRUE
# Mock      0.9997584541 0.0002415459
# Bystander 0.6937901499 0.3062098501
# Infected  0.0000000000 1.0000000000

table(sce$Library, sce$virus_pct > 50)
#           FALSE TRUE
# Mock       4140    0
# Bystander  2801    1
# Infected   1029 1358

sce$virus_pct[sce$Library == "Infected" ] %>% summary() 
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.3501 13.6756 62.2449 50.7004 82.7895 95.9093


#look at distributions on log10 scale----
#get half the minimum non-zero value

halfmin_virus <- min(sce$virus_pct[sce$virus_pct>0])
halfmin_virus
#0.002442659

sce$virus_pct_log10 <- log10(sce$virus_pct+halfmin_virus)
range(sce$virus_pct_log10)
#-2.612137  1.981872




x11(height = 10, width = 6)
layout(matrix(1:3,3,1))
hist(sce$virus_pct_log10[sce$Library=="Mock"], 1000, xlim = c(-3,2),
     main = "Mock virus percentage")
hist(sce$virus_pct_log10[sce$Library=="Bystander"],1000, xlim = c(-3,2),
     main = "Bystander virus percentage")
hist(sce$virus_pct_log10[sce$Library=="Infected"],1000, xlim = c(-3,2),
     main = "Infected virus percentage")

x11(height = 10, width = 6)
layout(matrix(1:3,3,1))
hist(sce$virus_pct_log10[sce$Library=="Mock"], 1000, xlim = c(-3,2),
     main = "Mock virus percentage", ylim = c(0,30))
hist(sce$virus_pct_log10[sce$Library=="Bystander"],1000, xlim = c(-3,2),
     main = "Bystander virus percentage", ylim = c(0,30))
hist(sce$virus_pct_log10[sce$Library=="Infected"],1000, xlim = c(-3,2),
     main = "Infected virus percentage", ylim = c(0,30))


#Find Infected threshold for virus_pct_log10  ----

#first calculate density and find local maxima
#taken from https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima

des.all <- density(sce$virus_pct_log10[sce$Library=="Infected"])
min.all <- des.all$x[which(diff(sign(diff(des.all$y)))==2)+1]
min.all
#0.4582323


x11(height = 10, width = 6)
jpeg("results/test_filter_and_normalize/2020-04-21-perth09_VirusPcts.jpeg",
     height = 10, width = 6, units = "in", res = 300, quality = 100)
layout(matrix(1:3,3,1))
par(mar = c(3,5,4,1))
hist(sce$virus_pct_log10[sce$Library=="Mock"], 1000, xlim = c(-3,2),
     main = "Mock", ylim = c(0,40),cex.lab = 2,cex.main = 2, cex.axis = 1.5,
     xlab = "", ylab = "cell number")
abline(v = min.all[1], col = 3, lty = 2, lwd = 2)
hist(sce$virus_pct_log10[sce$Library=="Bystander"],1000, xlim = c(-3,2),
     main = "Bystander", ylim = c(0,40),cex.lab = 2,cex.main = 2, cex.axis = 1.5,
     xlab = "", ylab = "cell number")
abline(v = min.all[1], col = 3, lty = 2, lwd = 2)
par(mar = c(6,5,4,1))
hist(sce$virus_pct_log10[sce$Library=="Infected"],1000, xlim = c(-3,2),
     main = "Infected", ylim = c(0,40),cex.lab = 2,cex.main = 2, cex.axis = 1.5,
     xlab = "log10(% of mRNA from flu)", ylab = "cell number")
abline(v = min.all[1], col = 3, lty = 2, lwd = 2)

dev.off()



#call infected any cell above the first minimum ----

sce$InfectedStatus <- "NotInfected"
sce$InfectedStatus[sce$virus_pct_log10 > min.all[1]] <- "Infected"

table(sce$Library, sce$InfectedStatus) 
#           Infected NotInfected
# Mock             0        4140
# Bystander        6        2796
# Infected      1940         447


table(sce$Library, sce$InfectedStatus)%>% prop.table(margin = 1)
#              Infected  NotInfected
# Mock      0.000000000 1.0000000000
# Bystander 0.002141328 0.997858672
# Infected  0.812735651 0.187264349



## Calculate individual virus gene percentages----

n.colData <- ncol(colData(sce)) 
virus_gene_counts <- tail(assay(sce, "counts"), 8) %>% t() %>% as.matrix()
colData(sce) <- cbind(colData(sce), virus_gene_counts)
names(colData(sce))[(n.colData+1):(n.colData+8)] <- paste0(names(colData(sce))[(n.colData+1):(n.colData+8)],"_counts")
n.colData <- ncol(colData(sce)) 
colData(sce) <- cbind(colData(sce), virus_gene_counts/sce$total_counts * 100)
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


vir_pcts <- colData(sce)[,grep("pct$", names(colData(sce)))[-1]]
vir_pcts_log10 <- colData(sce)[,grep("log10$", names(colData(sce)))[-1]]

x11(width = 6, height = 10)
layout(matrix(1:8, 4, 2, byrow = TRUE))
for (i in 1:8){
  hist(vir_pcts_log10[sce$Library == "Infected",i], 1000, xlim = c(-3,2), 
       main = names(vir_pcts_log10)[i])
}
#PB2, PB1, and PA have many zeros

#Try without zeros

x11(width = 6, height = 10)
layout(matrix(1:8, 4, 2, byrow = TRUE))
for (i in 1:8){
  hist(vir_pcts_log10[sce$Library == "Infected" & vir_pcts[,i] > 0,i], 
       1000, xlim = c(-3,2), 
       main = names(vir_pcts_log10)[i])
}


# Calculate minima for non-zero infected cells ----

min.indiv <- rep(NA, 8)
names(min.indiv) <- tail(rownames(sce), 8)

for (i in 1:8){
  #first calculate density and find local maxima
  des.1v <-density(vir_pcts_log10[sce$Library == "Infected" & vir_pcts[,i] > 0,i])
  min.indiv[i] <- des.1v$x[which(diff(sign(diff(des.1v$y)))==2)+1][1]
  jpeg(paste0("results/test_filter_and_normalize/2020-04-02-perth09_",names(min.indiv)[i],"_genes.jpeg"),
       height = 7, width = 7, units = "in", res = 300, quality = 100)
  layout(matrix(1:2,2,1))
  temp1 <- vir_pcts_log10[sce$Library == "Bystander",i]
  temp2 <- vir_pcts_log10[sce$Library == "Bystander" & vir_pcts[,i] > 0 & sce$InfectedStatus == "NotInfected",i]
  hist(temp1,1000, xlim = c(-3,2),
       main = paste("Bystander",names(vir_pcts)[i]), ylim = c(0,30))
  abline(v = quantile(temp2, 0.95), col = 2, lty = 2, lwd = 2)
  hist(vir_pcts_log10[sce$Library == "Infected",i],1000, xlim = c(-3,2),
       main = paste("Infected",names(vir_pcts)[i]), ylim = c(0,30))
  abline(v = min.indiv[i], col = 3, lty = 2, lwd = 2)
  dev.off()
}

#This worked well for everyone 

min.indiv
#        PB2        PB1         PA         HA         NP         NA          M 
# -1.4504983 -1.1141712 -1.5278023 -0.4161672 -0.1506843 -0.9012465  0.1165561 
#         NS 
# -0.2509913 


#Put on one plot

x11(width = 6, height = 10)
layout(matrix(1:8, 4, 2, byrow = TRUE))
for (i in 1:8){
  hist(vir_pcts_log10[sce$Library == "Infected" & vir_pcts[,i] > 0,i], 
       1000, xlim = c(-3,2), 
       main = names(vir_pcts_log10)[i])
  abline(v = min.indiv[i], col = 3, lty = 2, lwd = 2)
}
  



jpeg("results/test_filter_and_normalize/2020-04-17-perth09_8_genes.jpeg",
     height = 10, width = 6, units = "in", res = 300, quality = 100)
layout(matrix(1:8, 4, 2, byrow = TRUE))
for (i in 1:8){
  hist(vir_pcts_log10[sce$Library == "Infected" & vir_pcts[,i] > 0,i], 
       1000, xlim = c(-3,2), xlab = "log10(% of mRNA from each flu gene)", ylim = c(0,25), 
       main = gsub("_pct_log10","",names(vir_pcts_log10)[i]), ylab = "cell number")
  abline(v = min.indiv[i], col = 3, lty = 2, lwd = 2)
}
dev.off()


#call present/absent for each gene ----

sce$PB2_status <- ifelse(sce$PB2_pct_log10 > min.indiv["PB2"], "Present", "Absent")
table(sce$Library,sce$PB2_status)

sce$PB1_status <- ifelse(sce$PB1_pct_log10 > min.indiv["PB1"], "Present", "Absent")
table(sce$Library,sce$PB1_status)

sce$PA_status <- ifelse(sce$PA_pct_log10 > min.indiv["PA"], "Present", "Absent")
table(sce$Library,sce$PA_status)

sce$HA_status <- ifelse(sce$HA_pct_log10 > min.indiv["HA"], "Present", "Absent")
table(sce$Library,sce$HA_status)

sce$NP_status <- ifelse(sce$NP_pct_log10 > min.indiv["NP"], "Present", "Absent")
table(sce$Library,sce$NP_status)

sce$NA_status <- ifelse(sce$NA_pct_log10 > min.indiv["NA"], "Present", "Absent")
table(sce$Library,sce$NA_status)

sce$M_status <- ifelse(sce$M_pct_log10 > min.indiv["M"], "Present", "Absent")
table(sce$Library,sce$M_status)

sce$NS_status <- ifelse(sce$NS_pct_log10 > min.indiv["NS"], "Present", "Absent")
table(sce$Library,sce$NS_status)


#calculate the number viral genes present ----

temp <-  as.matrix(colData(sce)[,grep("_status$", names(colData(sce)))])
sce$NumPres <- rowSums(temp == "Present")
table(sce$Library, sce$NumPres)
#              0    1    2    3    4    5    6    7    8
# Mock      4140    0    0    0    0    0    0    0    0
# Bystander 2787    4    3    2    0    4    1    1    0
# Infected   417   20    9   20   40   90  245  529 1017


##Save rds file----

saveRDS(sce, file = "results/test_filter_and_normalize/2020-04-02-perth09_finalHostVirus.rds")
#sce <- readRDS("results/test_filter_and_normalize/2020-04-02-perth09_finalHostVirus.rds")

##Create metadata file ####
meta.1 <- colData(sce)
out.meta <- "results/test_filter_and_normalize/2020-04-02-perth09_metadata.tsv"
write.table(meta.1, file = out.meta, sep = "\t", row.names = FALSE)


##Subset to infected Library ----

sce.inf <- sce[,sce$Library == "Infected" & sce$InfectedStatus == "Infected"]

#Re filter genes to keep only ones present in at least 5 infected cells ----
#(this is same as what happens in Seurat testing)

dim(sce.inf)
#19235  1940
keep.feature <- nexprs(sce.inf, byrow=TRUE) >= 5
sce.inf <- sce.inf[keep.feature,]
remove(keep.feature)
dim(sce.inf)
#12877  1940


##Do normalization ----

clusters <- quickCluster(sce.inf, use.ranks=FALSE, BSPARAM=IrlbaParam())

##create size factors for normalizing within clusters
sce.inf <- computeSumFactors(sce.inf, min.mean=0.1, cluster=clusters)
sf <- sce.inf@int_colData@listData$size_factor
sce.inf <- computeSumFactors(sce.inf, min.mean=0.1,scaling = sf)
#remove(clusters)
#remove(sf)

#Remove viral genes (same as Seurat testing)

sce.inf <- sce.inf[!rownames(sce.inf) %in% c("PB2","PB1","PA","HA","NP","NA","M","NS"),]


##Save rds file for NBID testing ----

saveRDS(sce.inf, file = "results/test_filter_and_normalize/2020-04-18-perth09_InfectedOnly.rds")


## Import this into biocluster to do testing there



sessionInfo()
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 18363)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods  
# [9] base     
# 
# other attached packages:
#   [1] magrittr_1.5                BiocSingular_1.2.2          scran_1.14.6               
# [4] EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.10.2            AnnotationFilter_1.10.0    
# [7] GenomicFeatures_1.38.2      AnnotationDbi_1.48.0        scater_1.14.6              
# [10] ggplot2_3.3.0               DropletUtils_1.6.1          SingleCellExperiment_1.8.0 
# [13] SummarizedExperiment_1.16.1 DelayedArray_0.12.2         BiocParallel_1.20.1        
# [16] matrixStats_0.56.0          Biobase_2.46.0              GenomicRanges_1.38.0       
# [19] GenomeInfoDb_1.22.1         IRanges_2.20.2              S4Vectors_0.24.3           
# [22] BiocGenerics_0.32.0         simpleSingleCell_1.10.1    
# 
# loaded via a namespace (and not attached):
#   [1] ProtGenerics_1.18.0      bitops_1.0-6             bit64_0.9-7             
# [4] progress_1.2.2           httr_1.4.1               tools_3.6.3             
# [7] R6_2.4.1                 irlba_2.3.3              HDF5Array_1.14.3        
# [10] vipor_0.4.5              lazyeval_0.2.2           DBI_1.1.0               
# [13] colorspace_1.4-1         withr_2.1.2              tidyselect_1.0.0        
# [16] gridExtra_2.3            prettyunits_1.1.1        bit_1.1-15.2            
# [19] curl_4.3                 compiler_3.6.3           cli_2.0.2               
# [22] BiocNeighbors_1.4.2      rtracklayer_1.46.0       scales_1.1.0            
# [25] askpass_1.1              rappdirs_0.3.1           Rsamtools_2.2.3         
# [28] stringr_1.4.0            digest_0.6.25            R.utils_2.9.2           
# [31] XVector_0.26.0           pkgconfig_2.0.3          dbplyr_1.4.2            
# [34] limma_3.42.2             rlang_0.4.5              rstudioapi_0.11         
# [37] RSQLite_2.2.0            DelayedMatrixStats_1.8.0 dplyr_0.8.5             
# [40] R.oo_1.23.0              RCurl_1.98-1.1           GenomeInfoDbData_1.2.2  
# [43] Matrix_1.2-18            Rcpp_1.0.4.6             ggbeeswarm_0.6.0        
# [46] munsell_0.5.0            Rhdf5lib_1.8.0           fansi_0.4.1             
# [49] viridis_0.5.1            lifecycle_0.2.0          R.methodsS3_1.8.0       
# [52] stringi_1.4.6            edgeR_3.28.1             zlibbioc_1.32.0         
# [55] rhdf5_2.30.1             BiocFileCache_1.10.2     grid_3.6.3              
# [58] blob_1.2.1               dqrng_0.2.1              crayon_1.3.4            
# [61] lattice_0.20-38          Biostrings_2.54.0        hms_0.5.3               
# [64] locfit_1.5-9.4           pillar_1.4.3             igraph_1.2.5            
# [67] biomaRt_2.42.1           XML_3.99-0.3             glue_1.4.0              
# [70] vctrs_0.2.4              gtable_0.3.0             openssl_1.4.1           
# [73] purrr_0.3.3              assertthat_0.2.1         rsvd_1.0.3              
# [76] viridisLite_0.3.0        tibble_3.0.0             GenomicAlignments_1.22.1
# [79] beeswarm_0.2.3           memoise_1.1.0            statmod_1.4.34          
# [82] ellipsis_0.3.0  
 