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


#RESTART here 1 ----
sce <- readRDS("results/test_filter_and_normalize/2020-03-09_sce_postCellCycle.rds")

#### Add in library IDs ####

Lib.ID <- strsplit(sce$Barcode, "-") %>% sapply(function(x) x[2]) %>% factor
levels(Lib.ID) <- c("Bystander","Mock","Infected")
sce$Library <- factor(Lib.ID, levels = c("Mock","Bystander","Infected"))
table(sce$Library)
#Mock Bystander  Infected 
#4867      3758      4973


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

##Filter 3: Call doublets####
#Remove cells with more than twice the median total counts
#USING ONLY HOST COUNTS

#temp1 <- colSums(head(assay(sce, "normcounts"), -8))
temp1 <- colSums(head(assay(sce, "counts"), -8))
sce <- sce[,temp1 < 2*median(temp1)]
rm(temp1)
dim(sce)
#19779 12531
table(sce$Library)
#Mock Bystander  Infected 
#4314      3367      4850


##calculate final QC stats####
sce <- calculateQCMetrics(sce)  ###no MT control here, cells are infected


##add chromosome location info
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID,column="SEQNAME", keytype="GENEID")
rowData(sce)$chr <- location
remove(location)

#### Save copy for import into Seurat ####

out <- "results/test_filter_and_normalize/2020-03-30-cal07"
saveRDS(sce, file = paste0(out, "_sce_noNorm.rds"))

#RESTART here 2 ----
sce <- readRDS(paste0(out, "_sce_noNorm.rds"))



##DO NOT DO normalization####

#### Calculate host counts and viral counts ----

sce$host_counts <- colSums(head(assay(sce, "counts"), -8))
sce$virus_counts <- colSums(tail(assay(sce, "counts"), 8))
sce$virus_pct <- sce$virus_counts / sce$total_counts *100


table(sce$Library, sce$virus_pct >0) %>% prop.table(margin = 1)
#                  FALSE       TRUE
# Mock      0.93579045 0.06420955
# Bystander 0.94564895 0.05435105
# Infected  0.00000000 1.00000000


table(sce$Library, sce$virus_pct > 50)
#           FALSE TRUE
# Mock       4314    0
# Bystander  3367    0
# Infected   4638  212

sce$virus_pct[sce$Library == "Infected" ] %>% summary() 
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2424  5.0595 10.0196 13.8803 16.5209 92.1365




#look at distributions on log10 scale----
#get half the minimum non-zero value

halfmin_virus <- min(sce$virus_pct[sce$virus_pct>0])
halfmin_virus
#0.002834226
sce$virus_pct_log10 <- log10(sce$virus_pct+halfmin_virus)

range(sce$virus_pct_log10)
#-2.547566  1.964445


x11(height = 10, width = 6)
layout(matrix(1:3,3,1))
hist(sce$virus_pct_log10[sce$Library=="Mock"], 1000, xlim = c(-3,2),
     main = "Mock virus percentage", ylim = c(0,25))
hist(sce$virus_pct_log10[sce$Library=="Bystander"],1000, xlim = c(-3,2),
     main = "Bystander virus percentage", ylim = c(0,25))
hist(sce$virus_pct_log10[sce$Library=="Infected"],1000, xlim = c(-3,2),
     main = "Infected virus percentage", ylim = c(0,25))


#Find Infected threshold for virus_pct_log10  ----

#first calculate density and find local maxima
#taken from https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima

des.all <- density(sce$virus_pct_log10[sce$Library=="Infected"])
min.all <- des.all$x[which(diff(sign(diff(des.all$y)))==2)+1]
min.all
#-0.04843422  1.72627621

x11(height = 10, width = 6)
jpeg("results/test_filter_and_normalize/2020-04-02-cal07_VirusPcts.jpeg",
    height = 10, width = 6, units = "in", res = 300, quality = 100)
layout(matrix(1:3,3,1))
hist(sce$virus_pct_log10[sce$Library=="Mock"], 1000, xlim = c(-3,2),
     main = "Mock virus percentage", ylim = c(0,25),
     xlab = "log10(viral gene percentage)")
abline(v = quantile(sce$virus_pct_log10[sce$virus_pct>0 & sce$Library == "Mock"], 0.95), col = "blue", lty = 2, lwd = 2)
hist(sce$virus_pct_log10[sce$Library=="Bystander"],1000, xlim = c(-3,2),
     main = "Bystander virus percentage", ylim = c(0,25),
     xlab = "log10(viral gene percentage)")
abline(v = quantile(sce$virus_pct_log10[sce$virus_pct>0 & sce$Library == "Bystander"], 0.95), col = 2, lty = 2, lwd = 2)
hist(sce$virus_pct_log10[sce$Library=="Infected"],1000, xlim = c(-3,2),
     main = "Infected virus percentage", ylim = c(0,25),
     xlab = "log10(viral gene percentage)")
abline(v = min.all[1], col = 3, lty = 2, lwd = 2)

dev.off()


#call infected any cell above the first minimum ----

sce$InfectedStatus <- "NotInfected"
sce$InfectedStatus[sce$virus_pct_log10 > min.all[1]] <- "Infected"

table(sce$Library, sce$InfectedStatus) 
#           Infected NotInfected
# Mock             0        4314
# Bystander        2        3365
# Infected      4385         465


table(sce$Library, sce$InfectedStatus)%>% prop.table(margin = 1)
#                 Infected  NotInfected
# Mock      0.0000000000 1.0000000000
# Bystander 0.0005940006 0.9994059994
# Infected  0.9041237113 0.0958762887





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
  jpeg(paste0("results/test_filter_and_normalize/2020-04-02-cal07_",names(min.indiv)[i],"_genes.jpeg"),
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

#This worked well for everyone but PB2

min.indiv
#         PB2         PB1          PA          HA          NP          NA           M 
# -0.03778983 -1.50115874 -1.91880469 -0.72637811 -0.72424774 -1.22532475 -0.38626717 
#          NS 
# -0.68823119


#Put on one plot

x11(width = 6, height = 10)
layout(matrix(1:8, 4, 2, byrow = TRUE))
for (i in 1:8){
  hist(vir_pcts_log10[sce$Library == "Infected" & vir_pcts[,i] > 0,i], 
       1000, xlim = c(-3,2), 
       main = names(vir_pcts_log10)[i])
  abline(v = min.indiv[i], col = 3, lty = 2, lwd = 2)
}

# PB2 and PA don't look to good; these are low expression and so 
# set to -1.5

min.indiv[c(1,3)] <- -1.5


jpeg("results/test_filter_and_normalize/2020-04-02-cal07_8_genes.jpeg",
     height = 10, width = 6, units = "in", res = 300, quality = 100)
layout(matrix(1:8, 4, 2, byrow = TRUE))
for (i in 1:8){
  hist(vir_pcts_log10[sce$Library == "Infected" & vir_pcts[,i] > 0,i], 
       1000, xlim = c(-3,2), 
       main = names(vir_pcts_log10)[i])
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
# Mock      4314    0    0    0    0    0    0    0    0
# Bystander 3359    4    2    1    0    0    1    0    0
# Infected   406   55   34   82  106  204  593 1358 2012



##Save rds file----

saveRDS(sce, file = "results/test_filter_and_normalize/2020-04-02-cal07_finalHostVirus.rds")
#sce <- readRDS("results/test_filter_and_normalize/2020-04-02-cal07_finalHostVirus.rds")

##Create metadata file ####
meta.1 <- colData(sce)
out.meta <- "results/test_filter_and_normalize/2020-04-02-cal07_metadata.tsv"
write.table(meta.1, file = out.meta, sep = "\t", row.names = FALSE)


##Subset to infected Library ----

sce.inf <- sce[,sce$Library == "Infected"]


##Do normalization ----

clusters <- quickCluster(sce.inf, use.ranks=FALSE, BSPARAM=IrlbaParam())

##create size factors for normalizing within clusters
sce.inf <- computeSumFactors(sce.inf, min.mean=0.1, cluster=clusters)
sf <- sce.inf@int_colData@listData$size_factor
sce.inf <- computeSumFactors(sce.inf, min.mean=0.1,scaling = sf)
#remove(clusters)
#remove(sf)

##create logged & non-logged normalized values####
sce.inf <- logNormCounts(sce.inf,log = FALSE)
sce.inf <- logNormCounts(sce.inf,log = TRUE)


#Output raw and normalized counts for infected cells

normcounts <- as.matrix(assay(sce.inf,"normcounts"))
write.table(cbind(gene = rownames(normcounts), normcounts), 
            file = "results/test_filter_and_normalize/2020-04-03-cal07-Infected_Library_normcounts.txt",
            sep = "\t")

rawcounts <- as.matrix(assay(sce.inf,"counts"))
write.table(cbind(gene = rownames(rawcounts), rawcounts), 
            file = "results/test_filter_and_normalize/2020-04-03-cal07-Infected_Library_rawcounts.txt",
            sep = "\t")


#also subset to just InfectedStatus == Infected cells

normcounts <- normcounts[, sce.inf$InfectedStatus == "Infected"]
write.table(cbind(gene = rownames(normcounts), normcounts), 
            file = "results/test_filter_and_normalize/2020-04-03-cal07-Infected_Status_normcounts.txt",
            sep = "\t")

rawcounts <- rawcounts[, sce.inf$InfectedStatus == "Infected"]
write.table(cbind(gene = rownames(rawcounts), rawcounts), 
            file = "results/test_filter_and_normalize/2020-04-03-cal07-Infected_Status_rawcounts.txt",
            sep = "\t")


sessionInfo()
# sessionInfo()
# R version 3.6.2 (2019-12-12)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 18362)
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
#   [1] mixtools_1.2.0              magrittr_1.5                BiocSingular_1.2.1         
# [4] scran_1.14.5                EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.10.2           
# [7] AnnotationFilter_1.10.0     GenomicFeatures_1.38.0      AnnotationDbi_1.48.0       
# [10] scater_1.14.6               ggplot2_3.2.1               DropletUtils_1.6.1         
# [13] SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.1 DelayedArray_0.12.2        
# [16] BiocParallel_1.20.1         matrixStats_0.55.0          Biobase_2.46.0             
# [19] GenomicRanges_1.38.0        GenomeInfoDb_1.22.0         IRanges_2.20.2             
# [22] S4Vectors_0.24.2            BiocGenerics_0.32.0         simpleSingleCell_1.10.1    
# 
# loaded via a namespace (and not attached):
#   [1] segmented_1.1-0          ProtGenerics_1.18.0      bitops_1.0-6            
# [4] bit64_0.9-7              progress_1.2.2           httr_1.4.1              
# [7] tools_3.6.2              backports_1.1.5          R6_2.4.1                
# [10] irlba_2.3.3              HDF5Array_1.14.1         vipor_0.4.5             
# [13] DBI_1.1.0                lazyeval_0.2.2           colorspace_1.4-1        
# [16] withr_2.1.2              tidyselect_0.2.5         gridExtra_2.3           
# [19] prettyunits_1.1.0        bit_1.1-15.1             curl_4.3                
# [22] compiler_3.6.2           BiocNeighbors_1.4.1      rtracklayer_1.46.0      
# [25] scales_1.1.0             askpass_1.1              rappdirs_0.3.1          
# [28] stringr_1.4.0            digest_0.6.23            Rsamtools_2.2.1         
# [31] R.utils_2.9.2            XVector_0.26.0           pkgconfig_2.0.3         
# [34] dbplyr_1.4.2             limma_3.42.0             rlang_0.4.2             
# [37] rstudioapi_0.10          RSQLite_2.2.0            DelayedMatrixStats_1.8.0
# [40] dplyr_0.8.3              R.oo_1.23.0              RCurl_1.95-4.12         
# [43] GenomeInfoDbData_1.2.2   Matrix_1.2-18            Rcpp_1.0.3              
# [46] ggbeeswarm_0.6.0         munsell_0.5.0            Rhdf5lib_1.8.0          
# [49] viridis_0.5.1            lifecycle_0.1.0          R.methodsS3_1.7.1       
# [52] stringi_1.4.5            edgeR_3.28.0             MASS_7.3-51.4           
# [55] zlibbioc_1.32.0          rhdf5_2.30.1             BiocFileCache_1.10.2    
# [58] grid_3.6.2               blob_1.2.0               dqrng_0.2.1             
# [61] crayon_1.3.4             lattice_0.20-38          splines_3.6.2           
# [64] Biostrings_2.54.0        hms_0.5.3                locfit_1.5-9.1          
# [67] zeallot_0.1.0            pillar_1.4.3             igraph_1.2.4.2          
# [70] biomaRt_2.42.0           XML_3.98-1.20            glue_1.3.1              
# [73] vctrs_0.2.1              gtable_0.3.0             openssl_1.4.1           
# [76] purrr_0.3.3              kernlab_0.9-29           assertthat_0.2.1        
# [79] rsvd_1.0.2               survival_3.1-8           viridisLite_0.3.0       
# [82] tibble_2.1.3             GenomicAlignments_1.22.1 beeswarm_0.2.3          
# [85] memoise_1.1.0            statmod_1.4.33 
