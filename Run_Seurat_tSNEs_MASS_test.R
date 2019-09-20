#! /usr/bin/env Rscript

##Description: this R script imports a SingleCellExperiment object from SimpleSingleCell pipeline and runs Seurat tools
##             in order to perform DGE analysis and dimensional reduction analysis and produces DGE lists and tSNE
##             plots for several (default or from file) factors.


##libraries used/installed:
##Seurat_3.0.2
##simpleSingleCell_1.8.0
##DropletUtils_1.4.1
##ggplot2_3.1.1
##svglite_1.2.2
##dplyr_0.8.1
##magrittr_1.5

###defaults/settings
perp = 50          ##perplexity for tsne plot
do.plot1 = 0       ##generate first plot (top 10 bio-variance genes (non-virus))
do.plot2 = 0       ##generate second plot(s) (t-sne for various factors)
seed = 100         ##default seed

##specify which factors to plot in t-sne
factors2plot <- c("Library", "CellCycle", "TotalVirus", "TotalPB2", "TotalPB1", "TotalPA", "TotalHA", "TotalNP", "TotalNA", "TotalM", "TotalNS", "StatusInfected", "StatusPB2", "StatusPB1", "StatusPA", "StatusHA", "StatusNP", "StatusNA", "StatusM", "StatusNS", "MissingVirusGenes", "SingleMissingVirusGene", "NumMissingVirusGenes", "ClusterID", "Doublets")

##get command line
args <- commandArgs(TRUE)
infile <- args[1]          ###indir: input directory. this is a folder output from CellRanger containing the raw matrix (i.e. not the filtered matrix)
out <- args[2]            ###out: output base name (multiple files are output)
factorfile <- args[3]     ###factorfile: file containing single column of factors to produce tSNE plots for

##load initial SimpleSingleCell libraries
library(simpleSingleCell)
library(DropletUtils)
library(Seurat)
library(ggplot2)
library(svglite)
library(dplyr)
library(magrittr)

##import singlecellexperiment (sce) object
sce <- readRDS(infile)

##scan factor file (replace default factors)
if (file.exists(factorfile)){
  factors2plot <- scan(factorfile,what = "character")
}

#set default random seed
set.seed(seed)

##convert sce object to Seurat object (so)
so <- as.Seurat(sce, counts = "counts", data = "logcounts")

##scale data for PCA
all.genes <- rownames(so)
so <- ScaleData(so, features = all.genes)

##run PCA
so <- RunPCA(so, features = all.genes, seed.use = seed)

##perform clustering steps
so <- FindNeighbors(so, dims = 1:20)
so <- FindClusters(so, resolution = 0.5,random.seed = seed)

##run tSNE
so <- RunTSNE(object = so, dims = 1:20, perplexity = perp, seed.use = seed)

##subset SO object into infected only cells
so_infected <- subset(so, subset = StatusInfected == "Infected")

##redo scaling, PCA, and tSNE for infected cells
all.genes <- rownames(so_infected)
so_infected <- ScaleData(so_infected, features = all.genes)
so_infected <- RunPCA(so_infected, features = all.genes, seed.use = seed)
so_infected <- FindNeighbors(so_infected, dims = 1:20)
so_infected <- FindClusters(so_infected, resolution = 0.5,random.seed = seed)
so_infected <- RunTSNE(object = so_infected, dims = 1:20, perplexity = perp, seed.use = seed)

###run DGEs###
##run DGE between clusters for infected cells and print out a heatmap
so_infected.markers <- FindAllMarkers(so_infected, test.use = "MAST",  min.pct = 0.01, logfc.threshold = 0.25,random.seed = seed)
top10 <- so_infected.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_SeuratClusters.svg",sep = "")
write.table(so_infected.markers,file = out.file,sep = "\t")
plot <- DoHeatmap(so_infected, features = top10$gene) + NoLegend()
plot + theme(axis.text.y.left = element_text(size = rel(.5)))
out.file <- paste(out,"_Seurat_heatmap_InfectedCells_SeuratClusters.svg",sep = "")
ggsave(out.file, width = 8, height = 5)

##on Library: Infected vs DoubleNegative & Mock
Idents(so) <- "Library"
de <- FindMarkers(so, ident.1 = "Infected", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_AllCells_Library.txt",sep = "")
write.table(de,file = out.file,sep = "\t")
##on StatusInfected: Infected vs notInfected
Idents(so) <- "StatusInfected"
de <- FindMarkers(so, ident.1 = "Infected", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_AllCells_StatusInfected.txt",sep = "")
write.table(de,file = out.file,sep = "\t")

##individual missing virus genes
##for infected cells on PB2
Idents(so_infected) <- "StatusPB2"
de <- FindMarkers(so_infected, ident.1 = "Present", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_PB2.txt",sep = "")
write.table(de,file = out.file,sep = "\t")
##for infected cells on PB1
Idents(so_infected) <- "StatusPB1"
de <- FindMarkers(so_infected, ident.1 = "Present", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_PB1.txt",sep = "")
write.table(de,file = out.file,sep = "\t")
##for infected cells on PA
Idents(so_infected) <- "StatusPA"
de <- FindMarkers(so_infected, ident.1 = "Present", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_PA.txt",sep = "")
write.table(de,file = out.file,sep = "\t")
##for infected cells on HA
Idents(so_infected) <- "StatusHA"
de <- FindMarkers(so_infected, ident.1 = "Present", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_HA.txt",sep = "")
write.table(de,file = out.file,sep = "\t")
##for infected cells on NA
Idents(so_infected) <- "StatusNA"
de <- FindMarkers(so_infected, ident.1 = "Present", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_NA.txt",sep = "")
write.table(de,file = out.file,sep = "\t")
##for infected cells on NP
Idents(so_infected) <- "StatusNP"
de <- FindMarkers(so_infected, ident.1 = "Present", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_NP.txt",sep = "")
write.table(de,file = out.file,sep = "\t")
##for infected cells on M
Idents(so_infected) <- "StatusM"
de <- FindMarkers(so_infected, ident.1 = "Present", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_M.txt",sep = "")
write.table(de,file = out.file,sep = "\t")
##for infected cells on NS
Idents(so_infected) <- "StatusNS"
de <- FindMarkers(so_infected, ident.1 = "Present", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_NS.txt",sep = "")
write.table(de,file = out.file,sep = "\t")

###print out all cells tSNE plots###
##Seurat Clusters
plot <- DimPlot(so, reduction = "tsne")
plot + theme(title = element_text(size = rel(1.5)), plot.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.25)))
out.file <- paste(out,"_Seurat_tSNE_AllCells_SeuratClusters.svg",sep = "")
ggsave(out.file, width = 8, height = 5)
##library
plot <- DimPlot(so, reduction = "tsne",group.by = "Library")
plot + theme(title = element_text(size = rel(1.5)), plot.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.25)))
out.file <- paste(out,"_Seurat_tSNE_AllCells_Library.svg",sep = "")
ggsave(out.file, width = 8, height = 5)
##TotalVirus
plot <- FeaturePlot(so,reduction = "tsne", features = "TotalVirus")
plot + theme(title = element_text(size = rel(1.5)), plot.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.25)))
out.file <- paste(out,"_Seurat_tSNE_AllCells_TotalVirus.svg",sep = "")
ggsave(out.file, width = 8, height = 5)
##CellCycle
plot <- DimPlot(so, reduction = "tsne",group.by = "CellCycle")
plot + theme(title = element_text(size = rel(1.5)), plot.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.25)))
out.file <- paste(out,"_Seurat_tSNE_AllCells_CellCycle.svg",sep = "")
ggsave(out.file, width = 8, height = 5)


###print out infected cells tSNE plots###
##Seurat Clusters
plot <- DimPlot(so_infected, reduction = "tsne")
plot + theme(title = element_text(size = rel(1.5)), plot.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.25)))
out.file <- paste(out,"_Seurat_tSNE_InfectedCells_SeuratClusters.svg",sep = "")
ggsave(out.file, width = 8, height = 5)
##CellCycle
plot <- DimPlot(so_infected, reduction = "tsne",group.by = "CellCycle")
plot + theme(title = element_text(size = rel(1.5)), plot.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.25)))
out.file <- paste(out,"_Seurat_tSNE_InfectedCells_CellCycle.svg",sep = "")
ggsave(out.file, width = 8, height = 5)
##TotalVirus
plot <- FeaturePlot(so_infected,reduction = "tsne", features = "TotalVirus")
plot + theme(title = element_text(size = rel(1.5)), plot.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.25)))
out.file <- paste(out,"_Seurat_tSNE_InfectedCells_TotalVirus.svg",sep = "")
ggsave(out.file, width = 8, height = 5)

q()

