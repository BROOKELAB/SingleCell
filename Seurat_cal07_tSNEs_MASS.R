

###defaults/settings ----
perp = 50          ##perplexity for tsne plot
seed = 100         ##default seed

##specify which factors to plot in t-sne
factors2plot <- c("Library", "CellCycle", "TotalVirus", "TotalPB2", "TotalPB1", "TotalPA", "TotalHA", "TotalNP", "TotalNA", "TotalM", "TotalNS", "StatusInfected", "StatusPB2", "StatusPB1", "StatusPA", "StatusHA", "StatusNP", "StatusNA", "StatusM", "StatusNS", "MissingVirusGenes", "SingleMissingVirusGene", "NumMissingVirusGenes", "ClusterID", "Doublets")

##get command line
# args <- commandArgs(TRUE)
# infile <- args[1]          ###indir: input directory. this is a folder output from CellRanger containing the raw matrix (i.e. not the filtered matrix)
# out <- args[2]            ###out: output base name (multiple files are output)
# factorfile <- args[3]     ###factorfile: file containing single column of factors to produce tSNE plots for

##load initial SimpleSingleCell libraries
library(simpleSingleCell)
library(DropletUtils)
library(Seurat)
library(ggplot2)
library(svglite)
library(dplyr)
library(magrittr)

##import singlecellexperiment (sce) object ----
sce <- readRDS("results/test_filter_and_normalize/2020-04-02-cal07_finalHostVirus.rds")


sum(rowSums(assay(sce, "counts") >0 ) < 4)

#set default random seed
set.seed(seed)

##convert sce object to Seurat object (so)
#so <- as.Seurat(sce, counts = "counts")

#The above doesn't work because we didn't normalize the
#data. Instead, manuall make Seurat object:

so <- CreateSeuratObject(counts = as.matrix(assay(sce, "counts")), 
                         project = "cal07", meta.data = as.data.frame(colData(sce))[,-1])
so

head(so[[]])

rm(sce)
gc()

Idents(so) <- "Library"

x11(width = 12, height = 6)
VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA"), 
        ncol = 2, pt.size = 0.3)
#infected has lower numbers of gene detected
#The stark cutoff of nCount_RNA for Mock and Bystander but not infected came from
#the doublet filter on the non-viral genes only

x11()
FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.5 )


#### Normalization ####

# use the new sctransform method 
# https://satijalab.org/seurat/v3.1/sctransform_vignette.html

# Run first - it takes about 7 min

so
so <- SCTransform(so, return.only.var.genes = FALSE)
so

#A few more genes lost here because previous min_cell cutoff done before
#removing doublets


names(so)
# See new assay name
DefaultAssay(so)


# 
# ##scale data for PCA ----
# all.genes <- rownames(so)
# so <- ScaleData(so, features = all.genes)

##run PCA
all.genes <- rownames(so)
so <- RunPCA(so, features = all.genes, seed.use = seed)

x11()
ElbowPlot(so, 50)



##perform clustering steps
so <- FindNeighbors(so, dims = 1:40)
so <- FindClusters(so, resolution = 0.5,random.seed = seed)

##run tSNE ----
so <- RunTSNE(object = so, dims = 1:40, perplexity = perp, seed.use = seed)


#make figures for all cells ----

Idents(so) <- "seurat_clusters"
plot <- DimPlot(so, reduction = "tsne",pt.size = 1)
x11()
plot + theme(title = element_text(size = rel(1.5)), plot.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.25)))
out.file <- "results/Seurat_output/cal07_Seurat_tSNE_AllCells_SeuratClusters_2020-04-13.jpeg"
ggsave(out.file)

##library
#Add number of cells 
table(so$Library)
# Mock Bystander  Infected 
# 4314      3367      4850
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
plot <- DimPlot(so, reduction = "tsne",group.by = "Library",pt.size = 1)
x11(width = 9, height = 7)
plot + 
  scale_color_manual(values = ggplotColours(n = 3), labels = c("Mock (4314)", "Bystander (3367)", "Infected (4850)")) +
  theme(title = element_text(size = rel(1.5)), 
        plot.title = element_text(size = rel(2)), 
        legend.text = element_text(size = rel(1.75)),
        legend.position = c(0.65,0.95)) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
out.file <- "results/Seurat_output/cal07_Seurat_tSNE_AllCells_Library_2020-04-17.jpeg"
ggsave(out.file)

x11(width = 8, height = 18)
DimPlot(so, reduction = "tsne",group.by = "Library", split.by = "Library",pt.size = 1, ncol = 1) +
  theme(title = element_text(size = rel(1.5)), plot.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.25)))
out.file <- "results/Seurat_output/cal07_Seurat_tSNE_AllCells_Library_splitView_2020-04-13.jpeg"
ggsave(out.file)


plot <- DimPlot(so, reduction = "tsne",group.by = "InfectedStatus",pt.size = 1)
x11(width = 8, height = 7)
plot + theme(title = element_text(size = rel(1.5)), plot.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.25)))
out.file <- "results/Seurat_output/cal07_Seurat_tSNE_AllCells_InfectedStatus_2020-04-13.jpeg"
ggsave(out.file)


x11(width = 8, height = 18)
DimPlot(so, reduction = "tsne",group.by = "InfectedStatus", split.by = "Library",pt.size = 1, ncol = 1) +
  theme(title = element_text(size = rel(1.5)), plot.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.25)))
out.file <- "results/Seurat_output/cal07_Seurat_tSNE_AllCells_InfectedStatus_splitView_2020-04-13.jpeg"
ggsave(out.file)



## Try integration to remove batch effects 

#First, split the libraries into different objects
# see https://satijalab.org/seurat/v3.1/integration.html, SCTransform

# library.list <- SplitObject(so, split.by = "Library")
# for (i in 1:length(library.list)) {
#   library.list[[i]] <- SCTransform(library.list[[i]], verbose = FALSE, return.only.var.genes = FALSE)
# }
# 
# 
# library.features <- SelectIntegrationFeatures(object.list = library.list, nfeatures = 3000)
# 
# options(future.globals.maxSize = 3000 * 1024^2)
# library.list <- PrepSCTIntegration(object.list = library.list, anchor.features = library.features, 
#                                     verbose = FALSE)
# 
# library.anchors <- FindIntegrationAnchors(object.list = library.list, normalization.method = "SCT", 
#                                            anchor.features = library.features, verbose = FALSE)
# library.integrated <- IntegrateData(anchorset = library.anchors, normalization.method = "SCT", 
#                                      verbose = FALSE)
# 
# library.integrated <- RunPCA(library.integrated, verbose = FALSE)
# 
# library.integrated <- RunTSNE(object = library.integrated, dims = 1:40, perplexity = perp, seed.use = seed)
# 
# plot <- DimPlot(library.integrated, reduction = "tsne",group.by = "Library",pt.size = 1)
# x11(width = 8, height = 7)
# plot + theme(title = element_text(size = rel(1.5)), plot.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.25)))
# # out.file <- "results/Seurat_output/cal07_Seurat_tSNE_AllCells_Library.jpeg"
# # ggsave(out.file)
# 
# x11(width = 8, height = 18)
# DimPlot(library.integrated, reduction = "tsne",group.by = "Library", split.by = "Library",pt.size = 1, ncol = 1) +
#   theme(title = element_text(size = rel(1.5)), plot.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.25)))
# 
# #This removed all differences in Infected as well, so don't use



##subset SO object into infected only cells ----

so_infected <- subset(so, subset = Library == "Infected" & InfectedStatus == "Infected")


#Do SCTransform normalization
#NOTE: this will remove more genes because they aren't expressed in 
#just the infected cells

all.genes <- rownames(so_infected)
so_infected <- SCTransform(so_infected, return.only.var.genes = FALSE)
so_infected

#Do graph for Figure 2A ----

x11(5,5)

p <-ggplot(so_infected[[]], aes(x = virus_pct)) + 
  geom_histogram(fill = "blue") + 
  labs(x = "% of mRNA from flu", y = "Number of infected cells") +
  theme_classic(base_size = 20)
p
out.file <- "results/Seurat_output/cal07_Fig2A_2020-04-13.jpeg"
ggsave(out.file, width = 5, height = 5)

#get the bins and numbers

pg <- ggplot_build(p)
names(pg)
write.csv(pg$data, file = "results/Seurat_output/cal07_Fig2A_data_2020-04-17.csv",
          row.names = FALSE)

# vir_genes <- select(so_infected[[]], PB2_pct:NS_pct) %>% 
#   tibble::rownames_to_column("cell") %>%
#   tidyr::gather(PB2_pct:NS_pct, key = "gene", value = "pct") 
# 
# vir_genes$gene <- gsub("_pct", "", vir_genes$gene)
# vir_genes$gene <- factor(vir_genes$gene, level = unique(vir_genes$gene))
# 
# #add on half min non-zero value
# 
# vir_genes$pct <- vir_genes$pct + min(vir_genes$pct[vir_genes$pct > 0])
# 
# require(scales)
# x11(width = 9, height = 5)
# ggplot(vir_genes, aes(x=gene, y=pct, fill = gene)) + 
#   geom_violin(position = position_dodge(width = 0.9)) +
#   scale_y_continuous(trans='log10', limits = c(NA, 100), labels = comma) +
#   theme_classic(base_size = 20) +
#   labs(y = "Percentage of mRNA", x = "viral gene")  + NoLegend()
# #Not as good, but keep for now
# out.file <- "results/Seurat_output/cal07_Fig2B_2020-04-13.jpeg"
# ggsave(out.file, width = 9, height = 5)


#output normalized UMI values

temp <- as.matrix(so_infected[["SCT"]]@counts)

out.file <- "results/Seurat_output/cal07_NormalizedUMI_infected_2020-04-13.txt"
write.table(cbind(gene = rownames(temp), temp), file = out.file,
            row.names = FALSE, sep = "\t")

#add normalized viral gene values to metadata

temp2 <- t(tail(temp, n = 8))
colnames(temp2) <- paste0(colnames(temp2), "_normcounts")

so_infected@meta.data <- cbind(so_infected[[]], temp2)
names(so_infected[[]])

#don't write out meta.data now because the seurat cluster values are from
#all cells. Wait until the end.


#remove viral genes redo PCA, and tSNE for infected cells ----
newNum <- nrow(so_infected)-8
so_infected <- so_infected[1:newNum ,]
so_infected <- RunPCA(so_infected, features = all.genes, seed.use = seed)
so_infected <- FindNeighbors(so_infected, dims = 1:40)
so_infected <- FindClusters(so_infected, resolution = 0.5,random.seed = seed)
so_infected <- RunTSNE(object = so_infected, dims = 1:40, perplexity = perp, seed.use = seed)

Idents(so_infected) <- "seurat_clusters"
plot <- DimPlot(so_infected, reduction = "tsne",pt.size = 1)
x11(width = 8, height = 7)
plot + theme(title = element_text(size = rel(1.5)), plot.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.25)))
out.file <- "results/Seurat_output/cal07_Seurat_tSNE_Infected_SeuratClusters_2020-04-13.jpeg"
ggsave(out.file, width = 8, height = 7)



plot <- DimPlot(so_infected, reduction = "tsne",group.by = "CellCycle",pt.size = 1)
x11(width = 8, height = 7)
plot + theme(title = element_text(size = rel(1.5)), plot.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.25)))
out.file <- "results/Seurat_output/cal07_Seurat_tSNE_Infected_CellCycle_2020-04-13.jpeg"
ggsave(out.file)


plot <- FeaturePlot(so_infected,reduction = "tsne", features = "virus_pct",pt.size = 1)
x11(width = 8, height = 7)
plot + theme(title = element_text(size = rel(1.5)),  legend.text = element_text(size = rel(1.25)))
out.file <- "results/Seurat_output/cal07_Seurat_tSNE_Infected_ViralPct_2020-04-13.jpeg"
ggsave(out.file)




###run DGEs####
##run DGE between clusters for infected cells and print out a heatmap
so_infected.markers <- FindAllMarkers(so_infected, test.use = "MAST",  min.pct = 0.01, logfc.threshold = 0.25,random.seed = seed)
#select first 10 genes in each cluster (already in significance order)
top10 <- so_infected.markers %>% group_by(cluster) %>% filter(row_number() <= 10) %>% ungroup()
out.file <- "results/Seurat_output/cal07_Seurat_MAST_DGElist_InfectedCells_SeuratClusters_2020-04-13.txt"
write.table(so_infected.markers,file = out.file,sep = "\t")
plot <- DoHeatmap(so_infected, features = top10$gene) 
x11(width = 8, height = 5)
plot + theme(axis.text.y.left = element_text(size = rel(.5))) + NoLegend()
out.file <- "results/Seurat_output/cal07_Seurat_heatmap_InfectedCells_SeuratClusters_2020-04-17.jpeg"
ggsave(out.file, width = 8, height = 5)


# ##on Library: Infected vs DoubleNegative & Mock
# Idents(so) <- "Library"
# de <- FindMarkers(so, ident.1 = "Infected", test.use = "MAST",min.pct = 0.01, random.seed = seed)
# out.file <- paste(out,"_Seurat_MAST_DGElist_AllCells_Library.txt",sep = "")
# write.table(de,file = out.file,sep = "\t")
# ##on StatusInfected: Infected vs notInfected
# Idents(so) <- "StatusInfected"
# de <- FindMarkers(so, ident.1 = "Infected", test.use = "MAST",min.pct = 0.01, random.seed = seed)
# out.file <- paste(out,"_Seurat_MAST_DGElist_AllCells_StatusInfected.txt",sep = "")
# write.table(de,file = out.file,sep = "\t")


##individual missing virus genes
##for infected cells on PB2
out <- "results/Seurat_output/cal07"
Idents(so_infected) <- "PB2_status"
de <- FindMarkers(so_infected, ident.1 = "Present", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_PB2_2020-04-13.txt",sep = "")
write.table(de,file = out.file,sep = "\t")
##for infected cells on PB1
Idents(so_infected) <- "PB1_status"
de <- FindMarkers(so_infected, ident.1 = "Present", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_PB1_2020-04-13.txt",sep = "")
write.table(de,file = out.file,sep = "\t")
##for infected cells on PA
Idents(so_infected) <- "PA_status"
de <- FindMarkers(so_infected, ident.1 = "Present", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_PA_2020-04-13.txt",sep = "")
write.table(de,file = out.file,sep = "\t")
##for infected cells on HA
Idents(so_infected) <- "HA_status"
de <- FindMarkers(so_infected, ident.1 = "Present", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_HA_2020-04-13.txt",sep = "")
write.table(de,file = out.file,sep = "\t")
##for infected cells on NP
Idents(so_infected) <- "NP_status"
de <- FindMarkers(so_infected, ident.1 = "Present", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_NP.txt",sep = "")
write.table(de,file = out.file,sep = "\t")
##for infected cells on NA
Idents(so_infected) <- "NA_status"
de <- FindMarkers(so_infected, ident.1 = "Present", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_NA_2020-04-13.txt",sep = "")
write.table(de,file = out.file,sep = "\t")
##for infected cells on M
Idents(so_infected) <- "M_status"
de <- FindMarkers(so_infected, ident.1 = "Present", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_M_2020-04-13.txt",sep = "")
write.table(de,file = out.file,sep = "\t")
##for infected cells on NS
Idents(so_infected) <- "NS_status"
de <- FindMarkers(so_infected, ident.1 = "Present", test.use = "MAST",min.pct = 0.01, random.seed = seed)
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_NS_2020-04-13.txt",sep = "")
write.table(de,file = out.file,sep = "\t")


#output the metadata ----

out.meta <- "results/Seurat_output/cal07_MetaData_infected_2020-04-13.txt"
write.table(so_infected[[]], file = out.meta, sep = "\t", row.names = FALSE)





## Do tSNE of mock cells ----

so_mock <- subset(so, subset = Library == "Mock" & InfectedStatus == "NotInfected")


#Do SCTransform normalization
#NOTE: this will remove more genes because they aren't expressed in 
#just the mock cells

so_mock
so_mock <- SCTransform(so_mock, return.only.var.genes = FALSE)
so_mock


#remove viral genes redo PCA, and tSNE for mock cells ----
newNum <- nrow(so_mock)-8
so_mock <- so_mock[1:newNum ,]
so_mock <- RunPCA(so_mock, features = all.genes, seed.use = seed)
so_mock <- FindNeighbors(so_mock, dims = 1:40)
so_mock <- FindClusters(so_mock, resolution = 0.5,random.seed = seed)
so_mock <- RunTSNE(object = so_mock, dims = 1:40, perplexity = perp, seed.use = seed)

Idents(so_mock) <- "seurat_clusters"
plot <- DimPlot(so_mock, reduction = "tsne",pt.size = 1)
x11(width = 8, height = 7)
plot + theme(title = element_text(size = rel(1.5)), plot.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.25)))
out.file <- "results/Seurat_output/cal07_Seurat_tSNE_Mock_SeuratClusters_2020-04-13.jpeg"
ggsave(out.file, width = 8, height = 7)





## Output number of present and absent calls for each gene ----

table(so_infected$NumPres)
#   1    2    3    4    5    6    7    8 
#   2   26   81  105  204  593 1358 2012

names(so_infected[[]])

temp <- select(so_infected[[]], PB2_status:NS_status) %>% tidyr::gather(PB2_status:NS_status, key = "gene", value = "PresStatus")

table(temp$gene, temp$PresStatus)
#            Absent Present
# HA_status     273    4108
# M_status      413    3968
# NA_status     466    3915
# NP_status     232    4149
# NS_status     191    4190
# PA_status    1012    3369
# PB1_status    354    4027
# PB2_status   1210    3171

write.csv(table(so_infected$NumPres), 
          file = "results/Seurat_output/cal07_NumGenesPresPerCell_2020-04-17.csv")
write.csv(table(temp$gene, temp$PresStatus), 
          file = "results/Seurat_output/cal07_IndivGeneStatus_2020-04-17.csv")

save.image("results/Seurat_output/cal07.RData")
