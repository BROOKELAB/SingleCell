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
library(mixtools)


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


##Do normalization####
##generate quick clusters for normalization
clusters <- quickCluster(sce, use.ranks=FALSE, BSPARAM=IrlbaParam())

##create size factors for normalizing within clusters
sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters)
sce$sf <- sce@int_colData@listData$size_factor
sce <- computeSumFactors(sce, min.mean=0.1,scaling = sce$sf)
# remove(clusters)
# remove(sf)

range(sce$sf)
#0.0187505 2.7613344

table(sce$Library, clusters)


##create non-logged normalized values####
sce <- logNormCounts(sce,log = FALSE)
sce <- logNormCounts(sce,log = TRUE)

sce$diffcounts <- colSums(assay(sce, "normcounts")) - colSums(assay(sce, "counts"))

table(sce$diffcounts > 0)
#FALSE  TRUE 
# 4993  4336

table(sce$Library, sce$diffcounts > 0)
#           FALSE TRUE
# Mock       2773 1367
# Bystander  1910  892
# Infected    310 2077

tail(sce$diffcounts, 8)


#### Calculate total normcount, host normcounts and viral normcounts----

sce$total_normcounts <- colSums(assay(sce, "normcounts"))
sce$host_normcounts <- colSums(head(assay(sce, "normcounts"), -8))
sce$virus_normcounts <- colSums(tail(assay(sce, "normcounts"), 8))

sce$virus_pct <- sce$virus_normcounts / sce$total_normcounts *100

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


#Also check raw counts for sanity

sce$host_counts <- colSums(head(assay(sce, "counts"), -8))
sce$virus_counts <- colSums(tail(assay(sce, "counts"), 8))
sce$virus_pct_raw <- sce$virus_counts / sce$total_counts *100


#look at distributions on log10 scale----
#get half the minimum non-zero value

halfmin_virus <- min(sce$virus_pct[sce$virus_pct>0])
halfmin_virus
#0.002442659

sce$virus_pct_log10 <- log10(sce$virus_pct+halfmin_virus)
range(sce$virus_pct_log10)
#-2.612137  1.981872


# temp <- as.data.frame(colData(sce))
# 
# 
# 
# x11()
# ggplot(temp, aes(x=host_counts)) + 
#   geom_histogram(binwidth=1000) + 
#   facet_grid(Library ~ .)
# ggsave("results/test_filter_and_normalize/bioc_totalHostCounts.jpeg")
# #more infected cells have lower total host counts
# 
# #Do the libraries drastically differ in the number of non-zero genes?
# 
# x11()
# ggplot(temp, aes(x=total_features_by_counts)) + 
#   geom_histogram() + 
#   facet_grid(Library ~ .) 
# ggsave("results/test_filter_and_normalize/bioc_totalGenesDet.jpeg")
# #yes
# 
# x11()
# ggplot(temp, aes(x=host_normcounts)) + 
#   geom_histogram() + 
#   facet_grid(Library ~ .)
# ggsave("results/test_filter_and_normalize/bioc_totalHostNormCounts_scaling.jpeg")
# #why are regular gene counts shifted so far up?
# 
# x11()
# ggplot(temp, aes(x=virus_pct)) + 
#   geom_histogram(binwidth=10) + 
#   facet_grid(Library ~ .)
# 
# 
# x11()
# ggplot(temp, aes(x=virus_pct_log10)) + 
#   geom_histogram() + 
#   facet_grid(Library ~ .)
# 


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

#Calculate bimodal index and use + 3 SD of lower peak ----
#for virus_pct_log10 

#first calculate density and find local maxima
des.all <- density(sce$virus_pct_log10[sce$Library=="Infected"])

mu0.all <- des.all$x[which(diff(sign(diff(des.all$y)))==-2)+1]
mu0.all
# -0.166568  1.895880

#get inflection points
infl.all <- des.all$x[which(diff(sign(diff(des.all$y)))==2)+1]


#input this into the model to see the means:

mixmdl_virus = normalmixEM(sce$virus_pct_log10[sce$Library=="Infected"],
                           mu = mu0.all[1:2], lambda = c(0.1,0.9))
#NOTE: the above is not determinant and the exact values may
#depend on the number of iterations

x11()
plot(mixmdl_virus,which=2)
lines(des.all, lty=2, lwd=2)

mixmdl_virus$lambda
#[1] 0.5273708 0.4726292

#This isn't working well; try running with bottom 50% of cells

temp <- sce$virus_pct_log10[sce$Library=="Infected"] %>% sort() %>% head(n = round(sum(sce$Library=="Infected")/2))
#first calculate density and find local maxima
des.all <- density(temp)

mu0.all <- des.all$x[which(diff(sign(diff(des.all$y)))==-2)+1]

mu0.all
#-0.1502371  1.6193902


#input this into the model to see the means:

mixmdl_virus = normalmixEM(temp,
                           mu = mu0.all[1:2])
#NOTE: the above is not determinant and the exact values may
#depend on the number of iterations

x11()
plot(mixmdl_virus,which=2)
lines(des.all, lty=2, lwd=2)

mixmdl_virus$lambda
#0.369741 0.630259
mixmdl_virus$mu
#-0.1078127  1.3887098
mixmdl_virus$sigma
#0.1593826 0.3391466

i <- which.min(mixmdl_virus$mu)
thresh_infected <- mixmdl_virus$mu[i] - 3*mixmdl_virus$sigma[i]
thresh_high <- mixmdl_virus$mu[i] + 3*mixmdl_virus$sigma[i]

x11(height = 10, width = 6)
jpeg("results/test_filter_and_normalize/2020-03-30-Perth09_VirusPcts.jpeg",
     height = 10, width = 6, units = "in", res = 300, quality = 100)
layout(matrix(1:3,3,1))
hist(sce$virus_pct_log10[sce$Library=="Mock"], 1000, xlim = c(-3,2),
     main = "Mock virus percentage", ylim = c(0,30),
     xlab = "log10(viral gene percentage)")
abline(v = quantile(sce$virus_pct_log10[sce$virus_pct>0 & sce$Library == "Mock"], 0.95), col = "blue", lty = 2, lwd = 2)
hist(sce$virus_pct_log10[sce$Library=="Bystander"],1000, xlim = c(-3,2),
     main = "Bystander virus percentage", ylim = c(0,30),
     xlab = "log10(viral gene percentage)")
abline(v = quantile(sce$virus_pct_log10[sce$virus_pct>0 & sce$Library == "Bystander"], 0.95), col = 2, lty = 2, lwd = 2)
hist(sce$virus_pct_log10[sce$Library=="Infected"],1000, xlim = c(-3,2),
     main = "Infected virus percentage", ylim = c(0,30),
     xlab = "log10(viral gene percentage)")
abline(v = c(thresh_infected,thresh_high), col = 3, lty = 2, lwd = 2)
abline(v = infl.all[1], col = "purple", lty = 2, lwd = 2)
dev.off()


#call infected any cell that is > 3 SD----

sce$InfectedStatus <- "NotInfected"
sce$InfectedStatus[sce$virus_pct_log10> thresh_high] <- "Infected"

table(sce$Library, sce$InfectedStatus) 
#           Infected NotInfected
# Mock             0        4140
# Bystander        8        2794
# Infected      1943         444


table(sce$Library, sce$InfectedStatus)%>% prop.table(margin = 1)
#              Infected  NotInfected
# Mock      0.000000000 1.0000000000
# Bystander 0.002855103 0.997144897
# Infected  0.813992459 0.186007541



## Calculate individual virus gene percentages----

n.colData <- ncol(colData(sce)) 
virus_gene_counts <- tail(assay(sce, "counts"), 8) %>% t() %>% as.matrix()
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


#Try calc bimodal index of non-zero infected cells ----

for (i in 1:8){
  #first calculate density and find local maxima
  des.1v <-density(vir_pcts_log10[sce$Library == "Infected" & vir_pcts[,i] > 0,i])
  mu0.1v <- des.1v$x[which(diff(sign(diff(des.1v$y)))==-2)+1]
  min.pt <- des.1v$x[which(diff(sign(diff(des.1v$y)))==2)+1]
  #input this into the model to see the means; 
  #seed lambda and sigma with values from all 8 genes
  mixmdl.1v = normalmixEM(vir_pcts_log10[sce$Library == "Infected" & vir_pcts[,i] > 0,i],
                           mu = mu0.1v[1:2], lambda = c(0.1,0.9))
  j <- which.min(mixmdl.1v$mu)
  thresh_temp <- mixmdl.1v$mu[j] + 3*mixmdl.1v$sigma[j]
  
  x11(height = 7, width = 7)
  layout(matrix(1:2,2,1))
  temp1 <- vir_pcts_log10[sce$Library == "Bystander",i]
  temp2 <- vir_pcts_log10[sce$Library == "Bystander" & vir_pcts[,i] > 0 & sce$InfectedStatus == "NotInfected",i]
  hist(temp1,1000, xlim = c(-3,2),
       main = paste("Bystander",names(vir_pcts)[i]), ylim = c(0,30))
  abline(v = quantile(temp2, 0.95), col = 2, lty = 2, lwd = 2)
  hist(vir_pcts_log10[sce$Library == "Infected",i],1000, xlim = c(-3,2),
       main = paste("Infected",names(vir_pcts)[i]), ylim = c(0,30))
  abline(v = thresh_temp, col = 3, lty = 2, lwd = 2)
  abline(v = min.pt[1], col = 4, lty = 2, lwd = 2)
}
  


for (i in 1:8){
  #first calculate density and find local maxima
  des.1v <-density(vir_pcts_log10[sce$Library == "Infected" & vir_pcts[,i] > 0,i])
  mu0.1v <- des.1v$x[which(diff(sign(diff(des.1v$y)))==2)+1]
  thresh_temp <- mu0.1v[1]
  
  x11(height = 7, width = 7)
  layout(matrix(1:2,2,1))
  temp1 <- vir_pcts_log10[sce$Library == "Bystander",i]
  temp2 <- vir_pcts_log10[sce$Library == "Bystander" & vir_pcts[,i] > 0 & sce$InfectedStatus == "NotInfected",i]
  hist(temp1,1000, xlim = c(-3,2),
       main = paste("Bystander",names(vir_pcts)[i]), ylim = c(0,30))
  abline(v = quantile(temp2, 0.95), col = 4, lty = 2, lwd = 2)
  hist(vir_pcts_log10[sce$Library == "Infected",i],1000, xlim = c(-3,2),
       main = paste("Infected",names(vir_pcts)[i]), ylim = c(0,30))
  abline(v = thresh_temp, col = 3, lty = 2, lwd = 2)
  
}
  
  
  
  
# jpeg("results/test_filter_and_normalize/2020-03-30-Cal07_VirusPcts.jpeg",
#      height = 10, width = 6, units = "in", res = 300, quality = 100)
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
abline(v = c(thresh_infected,thresh_high), col = 3, lty = 2, lwd = 2)
dev.off()









#put on same y-axis

x11(width = 6, height = 10)
jpeg(paste0(out, "_viralgene_histograms.jpeg"), width = 6, height = 9, units = "in", quality=100, res = 300)
layout(matrix(1:8, 4, 2, byrow = TRUE))
for (i in grep("log10$", names(colData(sce)))[-1] ){
  hist(colData(sce)[sce$Library == "Infected",i], 1000, xlim = c(-3,2), 
       main = names(colData(sce))[i], ylim = c(0,15))
}
dev.off()


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



#Do bimodal indexes for all 8

mixmdl_PB2 = normalmixEM(sce$PB2_pct_log10[sce$Library == "Infected"])
mixmdl_PB1 = normalmixEM(sce$PB1_pct_log10[sce$Library == "Infected"])
mixmdl_PA = normalmixEM(sce$PA_pct_log10[sce$Library == "Infected"])
mixmdl_HA = normalmixEM(sce$HA_pct_log10[sce$Library == "Infected"])
mixmdl_NP = normalmixEM(sce$NP_pct_log10[sce$Library == "Infected"])
mixmdl_NA = normalmixEM(sce$NA_pct_log10[sce$Library == "Infected"])
mixmdl_M = normalmixEM(sce$M_pct_log10[sce$Library == "Infected"])
mixmdl_NS = normalmixEM(sce$NS_pct_log10[sce$Library == "Infected"])



x11(6,9)
jpeg(paste0(out, "_bimodality_curves.jpeg"), width = 6, height = 9, units = "in", quality=100, res = 300)
layout(matrix(1:8,4,2, byrow = TRUE))
plot(mixmdl_PB2,which=2, sub = "PB2 bimodality")
lines(density(sce$PB2_pct_log10[sce$Library == "Infected"]), lty=2, lwd=2)
plot(mixmdl_PB1,which=2, sub = "PB1 bimodality")
lines(density(sce$PB1_pct_log10[sce$Library == "Infected"]), lty=2, lwd=2)
plot(mixmdl_PA,which=2, sub = "PA bimodality")
lines(density(sce$PA_pct_log10[sce$Library == "Infected"]), lty=2, lwd=2)
plot(mixmdl_HA,which=2, sub = "HA bimodality")
lines(density(sce$HA_pct_log10[sce$Library == "Infected"]), lty=2, lwd=2)
plot(mixmdl_NP,which=2, sub = "NP bimodality")
lines(density(sce$NP_pct_log10[sce$Library == "Infected"]), lty=2, lwd=2)
plot(mixmdl_NA,which=2, sub = "NA bimodality")
lines(density(sce$NA_pct_log10[sce$Library == "Infected"]), lty=2, lwd=2)
plot(mixmdl_M,which=2, sub = "M bimodality")
lines(density(sce$M_pct_log10[sce$Library == "Infected"]), lty=2, lwd=2)
plot(mixmdl_NS,which=2, sub = "NS bimodality")
lines(density(sce$NS_pct_log10[sce$Library == "Infected"]), lty=2, lwd=2)
dev.off()






#subset to DoubleNegative (Bystander)----

sce.by <- sce[, sce$Library %in% "DoubleNegative"]
sce.by$NS_proportion <- sce.by$NS / sce.by$total_count


#For PB2, PB1, PA and NA, get 95% non-zero percentages of Bystander. But check for
#outliers that might be infected

#Do PB2 ----

by_PB2 <- sce.by$PB2 / sce.by$total_count
sum(by_PB2 > 0)
#18
which.max(by_PB2)
#612
x11()
hist(by_PB2, 1000, ylim = c(0, 10))

by_PB2_non0 <- by_PB2[by_PB2>0]
by_PB2_non0 <- sort(by_PB2_non0)

#Compare with infected PB2 to see if any bystander obviously infected

#log10 scale
halfmin_PB2 <- min(infectInfo.sm$PB2_proportion[infectInfo.sm$PB2_proportion>0])/2
halfmin_PB2

x11()
jpeg(paste0(out, "_PB2_thresholds.jpeg"), width = 7, height = 7, units = "in", res = 300, quality = 100)
layout(matrix(1:2,2,1))
temp <- hist(log10(infectInfo.sm$PB2_proportion + halfmin_PB2), 1000, 
             main = "Infected cells, PB2 proportion")
abline(v = log10(quantile(by_PB2_non0, 0.95)+halfmin_PB2), col = 2, lty = 2, lwd = 2)
abline(v = log10(quantile(head(by_PB2_non0, -1), 0.95)+halfmin_PB2), col = 3, lty = 2, lwd = 2)
hist(log10(by_PB2_non0 + halfmin_PB2), temp$breaks,
     main = paste(length(by_PB2_non0),"non-zero bystander cells"))
#one obviously infected
abline(v = log10(quantile(by_PB2_non0, 0.95)+halfmin_PB2), col = 2, lty = 2, lwd = 2)
abline(v = log10(quantile(head(by_PB2_non0, -1), 0.95)+halfmin_PB2), col = 3, lty = 2, lwd = 2)
dev.off()

#compare to previous smalled proportion that was present

quantile(by_PB2_non0, 0.95)
#        95% 
#0.001466067

quantile(head(by_PB2_non0, -1), 0.95)
#0.0002989828

min(infectInfo.sm$PB2_proportion[infectInfo.sm$PB2_Status %in% "Present"])
#0.000299

#new 95% without outlier is very close



#DO for PB1 ----

by_PB1 <- sce.by$PB1 / sce.by$total_count
sum(by_PB1 > 0)
#33
which.max(by_PB1)
#612
x11()
hist(by_PB1, 1000, ylim = c(0, 10))
#one outlier

by_PB1_non0 <- by_PB1[by_PB1>0]
by_PB1_non0 <- sort(by_PB1_non0)

#Compare with infected PB1 to see if any bystander obviously infected

#log10 scale
halfmin_PB1 <- min(infectInfo.sm$PB1_proportion[infectInfo.sm$PB1_proportion>0])/2
halfmin_PB1

x11()
jpeg(paste0(out, "_PB1_thresholds.jpeg"), width = 7, height = 7, units = "in", res = 300, quality = 100)
layout(matrix(1:2,2,1))
temp <- hist(log10(infectInfo.sm$PB1_proportion + halfmin_PB1), 1000, 
             main = "Infected cells, PB1 proportion")
abline(v = log10(quantile(by_PB1_non0, 0.95)+halfmin_PB1), col = 2, lty = 2, lwd = 2)
abline(v = log10(quantile(head(by_PB1_non0, -1), 0.95)+halfmin_PB1), col = 3, lty = 2, lwd = 2)
hist(log10(by_PB1_non0 + halfmin_PB1), temp$breaks,
     main = paste(length(by_PB1_non0),"non-zero bystander cells"))
#one obviously infected
abline(v = log10(quantile(by_PB1_non0, 0.95)+halfmin_PB1), col = 2, lty = 2, lwd = 2)
abline(v = log10(quantile(head(by_PB1_non0, -1), 0.95)+halfmin_PB1), col = 3, lty = 2, lwd = 2)
dev.off()

#compare to previous smalled proportion that was present

quantile(by_PB1_non0, 0.95)
#        95% 
#0.0002869663

quantile(head(by_PB1_non0, -1), 0.95)
#0.0002055767

min(infectInfo.sm$PB1_proportion[infectInfo.sm$PB1_Status %in% "Present"])
#0.000206

#new 95% without outlier is very close


#DO for PA ----

by_PA <- sce.by$PA / sce.by$total_count
sum(by_PA > 0)
#12
which.max(by_PA)
#612
x11()
hist(by_PA, 1000, ylim = c(0, 10))
#one outlier

by_PA_non0 <- by_PA[by_PA>0]
by_PA_non0 <- sort(by_PA_non0)

#Compare with infected PA to see if any bystander obviously infected

#log10 scale
halfmin_PA <- min(infectInfo.sm$PA_proportion[infectInfo.sm$PA_proportion>0])/2
halfmin_PA

x11()
jpeg(paste0(out, "_PA_thresholds.jpeg"), width = 7, height = 7, units = "in", res = 300, quality = 100)
layout(matrix(1:2,2,1))
temp <- hist(log10(infectInfo.sm$PA_proportion + halfmin_PA), 1000, 
             main = "Infected cells, PA proportion")
abline(v = log10(quantile(by_PA_non0, 0.95)+halfmin_PA), col = 2, lty = 2, lwd = 2)
abline(v = log10(quantile(head(by_PA_non0, -1), 0.95)+halfmin_PA), col = 3, lty = 2, lwd = 2)
hist(log10(by_PA_non0 + halfmin_PA), temp$breaks,
     main = paste(length(by_PA_non0),"non-zero bystander cells"))
#one obviously infected
abline(v = log10(quantile(by_PA_non0, 0.95)+halfmin_PA), col = 2, lty = 2, lwd = 2)
abline(v = log10(quantile(head(by_PA_non0, -1), 0.95)+halfmin_PA), col = 3, lty = 2, lwd = 2)
dev.off()

#compare to previous smalled proportion that was present

quantile(by_PA_non0, 0.95)
#        95% 
#0.002308459

quantile(head(by_PA_non0, -1), 0.95)
#0.0001367739

min(infectInfo.sm$PA_proportion[infectInfo.sm$PA_Status %in% "Present"])
#0.000137

#new 95% without outlier is very close


#DO for NA ----

by_NA <- sce.by$NA. / sce.by$total_count
sum(by_NA > 0)
#16
which.max(by_NA)
#563
x11()
hist(by_NA, 1000, ylim = c(0, 10))
#one outlier

by_NA_non0 <- by_NA[by_NA>0]
by_NA_non0 <- sort(by_NA_non0)

#Compare with infected NA to see if any bystander obviously infected

#log10 scale
halfmin_NA <- min(infectInfo.sm$NA_proportion[infectInfo.sm$NA_proportion>0])/2
halfmin_NA

x11()
jpeg(paste0(out, "_NA_thresholds.jpeg"), width = 7, height = 7, units = "in", res = 300, quality = 100)
layout(matrix(1:2,2,1))
temp <- hist(log10(infectInfo.sm$NA_proportion + halfmin_NA), 1000, 
             main = "Infected cells, NA proportion")
abline(v = log10(quantile(by_NA_non0, 0.95)+halfmin_NA), col = 2, lty = 2, lwd = 2)
abline(v = log10(quantile(head(by_NA_non0, -1), 0.95)+halfmin_NA), col = 3, lty = 2, lwd = 2)
hist(log10(by_NA_non0 + halfmin_NA), temp$breaks,
     main = paste(length(by_NA_non0),"non-zero bystander cells"))
#one obviously infected
abline(v = log10(quantile(by_NA_non0, 0.95)+halfmin_NA), col = 2, lty = 2, lwd = 2)
abline(v = log10(quantile(head(by_NA_non0, -1), 0.95)+halfmin_NA), col = 3, lty = 2, lwd = 2)
dev.off()

#compare to previous smalled proportion that was present

quantile(by_NA_non0, 0.95)
#        95% 
#0.001505873

quantile(head(by_NA_non0, -1), 0.95)
#0.0004488048

min(infectInfo.sm$NA_proportion[infectInfo.sm$NA_Status %in% "Present"])
#0.00045

#new 95% without outlier is very close


# For HA, NP, M and NS do something different----
#Calculate bimodal indexes and use + 3 SD of lower peak

library(mixtools)


#Start with HA----

by_HA <- sce.by$HA / sce.by$total_count
sum(by_HA > 0)
#45
which.max(by_HA)
#612
x11()
hist(by_HA, 1000, ylim = c(0, 10))
#one outlier

by_HA_non0 <- by_HA[by_HA>0]
by_HA_non0 <- sort(by_HA_non0)


#log10 scale
halfmin_HA <- min(infectInfo.sm$HA_proportion[infectInfo.sm$HA_proportion>0])/2


#need genes in rows and samples in columns

mixmdl_HA = normalmixEM(log10(infectInfo$HA_proportion+halfmin_HA))
x11()
plot(mixmdl_HA,which=2)
lines(density(log10(infectInfo$HA_proportion+halfmin_HA)), lty=2, lwd=2)

mixmdl_HA$mu
#-3.098479 -1.757742
mixmdl_HA$sigma
#0.1379327 0.4125467


#x11()
jpeg(paste0(out, "_HA_thresholds.jpeg"), width = 7, height = 7, units = "in", quality=100, res = 300)
layout(matrix(1:2,2,1))
temp <- hist(log10(infectInfo$HA_proportion+halfmin_HA), 1000,  
             main = "Infected cells, HA proportion")
abline(v = mixmdl_HA$mu[1] + 3*mixmdl_HA$sigma[1], col = 3, lty = 2, lwd = 2)
abline(v = log10(quantile(head(by_HA_non0, -1), 0.95)+halfmin_HA), col = 2, lty = 2, lwd = 2)
hist(log10(by_HA_non0 + halfmin_HA), temp$breaks,
     main = paste(length(by_HA_non0),"non-zero bystander cells"))
abline(v = mixmdl_HA$mu[1] + 3*mixmdl_HA$sigma[1], col = 3, lty = 2, lwd = 2)
abline(v = log10(quantile(head(by_HA_non0, -1), 0.95)+halfmin_HA), col = 2, lty = 2, lwd = 2)
dev.off()

#compare to previous smalled proportion that was present

quantile(by_HA_non0, 0.95)
#        95% 
#0.0008589049

quantile(head(by_HA_non0, -1), 0.95)
#0.0006606463

min(infectInfo.sm$HA_proportion[infectInfo.sm$HA_Status %in% "Present"])
#0.000661


#Do NP----

by_NP <- sce.by$NP / sce.by$total_count
sum(by_NP > 0)
#66
which.max(by_NP)
#563
x11()
hist(by_NP, 1000, ylim = c(0, 10))
#maybe two outlier?
which.max(by_NP[-563])
#611  - same as 612-1!

by_NP_non0 <- by_NP[by_NP>0]
by_NP_non0 <- sort(by_NP_non0)


#log10 scale
halfmin_NP <- min(infectInfo.sm$NP_proportion[infectInfo.sm$NP_proportion>0])/2


#need genes in rows and samples in columns

mixmdl_NP = normalmixEM(log10(infectInfo$NP_proportion+halfmin_NP))
x11()
plot(mixmdl_NP,which=2)
lines(density(log10(infectInfo$NP_proportion+halfmin_NP)), lty=2, lwd=2)

mixmdl_NP$mu
#-3.041578 -1.849426
mixmdl_NP$sigma
#0.1244391 0.4075653


#x11()
jpeg(paste0(out, "_NP_thresholds.jpeg"), width = 7, height = 7, units = "in", quality=100, res = 300)
layout(matrix(1:2,2,1))
temp <- hist(log10(infectInfo$NP_proportion+halfmin_NP), 1000,  
             main = "Infected cells, NP proportion", xlim = c(-4, -0.4))
abline(v = mixmdl_NP$mu[1] + 3*mixmdl_NP$sigma[1], col = 3, lty = 2, lwd = 2)
abline(v = log10(quantile(head(by_NP_non0, -1), 0.95)+halfmin_NP), col = 2, lty = 2, lwd = 2)
hist(log10(by_NP_non0 + halfmin_NP), 1000, xlim = c(-4, -0.4),
     main = paste(length(by_NP_non0),"non-zero bystander cells"))
abline(v = mixmdl_NP$mu[1] + 3*mixmdl_NP$sigma[1], col = 3, lty = 2, lwd = 2)
abline(v = log10(quantile(head(by_NP_non0, -1), 0.95)+halfmin_NP), col = 2, lty = 2, lwd = 2)
dev.off()

#compare to previous smalled proportion that was present

quantile(by_NP_non0, 0.95)
#        95% 
#0.0008485093

quantile(head(by_NP_non0, -1), 0.95)
#0.0006701869

min(infectInfo.sm$NP_proportion[infectInfo.sm$NP_Status %in% "Present"])
#0.000671



#Do NS----

by_NS <- sce.by$NS / sce.by$total_count
sum(by_NS > 0)
#80
which.max(by_NS)
#103
x11()
hist(by_NS, 1000, ylim = c(0, 10))
#maybe two outlier?
which.max(by_NS[-103])
#562

by_NS_non0 <- by_NS[by_NS>0]
by_NS_non0 <- sort(by_NS_non0)


#log10 scale
halfmin_NS <- min(infectInfo.sm$NS_proportion[infectInfo.sm$NS_proportion>0])/2


#need genes in rows and samples in columns

mixmdl_NS = normalmixEM(log10(infectInfo$NS_proportion+halfmin_NS))
x11()
plot(mixmdl_NS,which=2)
lines(density(log10(infectInfo$NS_proportion+halfmin_NS)), lty=2, lwd=2)

mixmdl_NS$mu
#-3.079796 -1.524843
mixmdl_NS$sigma
#0.1352574 0.3954511


#x11()
jpeg(paste0(out, "_NS_thresholds.jpeg"), width = 7, height = 7, units = "in", quality=100, res = 300)
layout(matrix(1:2,2,1))
temp <- hist(log10(infectInfo$NS_proportion+halfmin_NS), 1000,  
             main = "Infected cells, NS proportion", xlim = c(-4, -0.4))
abline(v = mixmdl_NS$mu[1] + 3*mixmdl_NS$sigma[1], col = 3, lty = 2, lwd = 2)
abline(v = log10(quantile(head(by_NS_non0, -1), 0.95)+halfmin_NS), col = 2, lty = 2, lwd = 2)
hist(log10(by_NS_non0 + halfmin_NS), 1000, xlim = c(-4, -0.4),
     main = paste(length(by_NS_non0),"non-zero bystander cells"))
abline(v = mixmdl_NS$mu[1] + 3*mixmdl_NS$sigma[1], col = 3, lty = 2, lwd = 2)
abline(v = log10(quantile(head(by_NS_non0, -1), 0.95)+halfmin_NS), col = 2, lty = 2, lwd = 2)
dev.off()

#compare to previous smalled proportion that was present

quantile(by_NS_non0, 0.95)
#        95% 
#0.001725732

quantile(head(by_NS_non0, -1), 0.95)
#0.0009229701

min(infectInfo.sm$NS_proportion[infectInfo.sm$NS_Status %in% "Present"])
#0.000923



#Do M----

by_M <- sce.by$M / sce.by$total_count
sum(by_M > 0)
#67
which.max(by_M)
#612
x11()
hist(by_M, 1000, ylim = c(0, 10))
#one outlier

by_M_non0 <- by_M[by_M>0]
by_M_non0 <- sort(by_M_non0)


#log10 scale
halfmin_M <- min(infectInfo.sm$M_proportion[infectInfo.sm$M_proportion>0])/2


#need genes in rows and samples in columns

mixmdl_M = normalmixEM(log10(infectInfo$M_proportion+halfmin_M))
x11()
plot(mixmdl_M,which=2)
lines(density(log10(infectInfo$M_proportion+halfmin_M)), lty=2, lwd=2)

mixmdl_M$mu
#-2.622773 -1.620764
mixmdl_M$sigma
#0.07911423 0.41559711


#x11()
jpeg(paste0(out, "_M_thresholds.jpeg"), width = 7, height = 7, units = "in", quality=100, res = 300)
layout(matrix(1:2,2,1))
temp <- hist(log10(infectInfo$M_proportion+halfmin_M), 1000,  
             main = "Infected cells, M proportion", xlim = c(-3.5, -0.2))
abline(v = mixmdl_M$mu[1] + 3*mixmdl_M$sigma[1], col = 3, lty = 2, lwd = 2)
abline(v = log10(quantile(head(by_M_non0, -1), 0.95)+halfmin_M), col = 2, lty = 2, lwd = 2)
hist(log10(by_M_non0 + halfmin_M), 1000, xlim = c(-3.5, -0.2),
     main = paste(length(by_M_non0),"non-zero bystander cells"))
abline(v = mixmdl_M$mu[1] + 3*mixmdl_M$sigma[1], col = 3, lty = 2, lwd = 2)
abline(v = log10(quantile(head(by_M_non0, -1), 0.95)+halfmin_M), col = 2, lty = 2, lwd = 2)
dev.off()

#compare to previous smalled proportion that was present

quantile(by_M_non0, 0.95)
#        95% 
#0.0007927077

quantile(head(by_M_non0, -1), 0.95)
#0.0005660847

min(infectInfo.sm$M_proportion[infectInfo.sm$M_Status %in% "Present"])
#0.000845
log10(0.000845++ halfmin_M)
#-2.897052


#Do bimodal indexes for all 8

mixmdl_PB2 = normalmixEM(log10(infectInfo$PB2_proportion+halfmin_PB2))
mixmdl_PB1 = normalmixEM(log10(infectInfo$PB1_proportion+halfmin_PB1))
mixmdl_PA = normalmixEM(log10(infectInfo$PA_proportion+halfmin_PA))
mixmdl_NA = normalmixEM(log10(infectInfo$NA_proportion+halfmin_NA))




x11(6,9)
jpeg(paste0(out, "_bimodality_curves.jpeg"), width = 6, height = 9, units = "in", quality=100, res = 300)
layout(matrix(1:8,4,2, byrow = TRUE))
plot(mixmdl_PB2,which=2, sub = "PB2 bimodality")
lines(density(log10(infectInfo$PB2_proportion+halfmin_PB2)), lty=2, lwd=2)
plot(mixmdl_PB1,which=2, sub = "PB1 bimodality")
lines(density(log10(infectInfo$PB1_proportion+halfmin_PB1)), lty=2, lwd=2)
plot(mixmdl_PA,which=2, sub = "PA bimodality")
lines(density(log10(infectInfo$PA_proportion+halfmin_PA)), lty=2, lwd=2)
plot(mixmdl_HA,which=2, sub = "HA bimodality")
lines(density(log10(infectInfo$HA_proportion+halfmin_HA)), lty=2, lwd=2)
plot(mixmdl_NP,which=2, sub = "NP bimodality")
lines(density(log10(infectInfo$NP_proportion+halfmin_NP)), lty=2, lwd=2)
plot(mixmdl_NA,which=2, sub = "NA bimodality")
lines(density(log10(infectInfo$NA_proportion+halfmin_NA)), lty=2, lwd=2)
plot(mixmdl_M,which=2, sub = "M bimodality")
lines(density(log10(infectInfo$M_proportion+halfmin_M)), lty=2, lwd=2)
plot(mixmdl_NS,which=2, sub = "NS bimodality")
lines(density(log10(infectInfo$NS_proportion+halfmin_NS)), lty=2, lwd=2)
dev.off()




##Save rds file----

saveRDS(sce, file = paste0(out, "_sce_finalHostVirus.rds"))
#sce <- readRDS(paste0(out, "_sce_finalHostVirus.rds"))

##Create first metadata file ####
meta.1 <- colData(sce)[,c("Barcode","Library","total_features_by_counts","total_counts","CellCycle")]
out.meta <- paste(out,"_1st_metadata.tsv",sep = "")
write.table(meta.1, file = out.meta, sep = "\t")




q()
