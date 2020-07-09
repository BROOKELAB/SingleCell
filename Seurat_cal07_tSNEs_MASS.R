

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


#set default random seed
set.seed(seed)

#Make Seurat object:

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
plot + theme(title = element_text(size = rel(1.5)),  legend.text = element_text(size = rel(1.25))) + labs(title = "")
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
plot + theme(axis.text.y.left = element_text(size = rel(.5))) + guides(color = FALSE)
out.file <- "results/Seurat_output/cal07_Seurat_heatmap_InfectedCells_SeuratClusters_2020-04-21.jpeg"
ggsave(out.file, width = 8, height = 5)



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
out.file <- paste(out,"_Seurat_MAST_DGElist_InfectedCells_NP_2020-04-13.txt",sep = "")
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


#Jiayi also wants me to add Seurat cluster number and NS status as separate
#rows in cal07_NormalizedUMI_infected_2020-04-13.txt

temp <- read.delim("results/Seurat_output/cal07_NormalizedUMI_infected_2020-04-13.txt")

#dashes in barcodes got changed to dots. change back:

names(temp) <- gsub(".","-", names(temp), fixed = TRUE)

temp2 <- t(so_infected[[]][,c("seurat_clusters","NS_status")]) %>% as.data.frame() %>%
          tibble::rownames_to_column(var = "gene")

all.equal(names(temp), names(temp2))


out.file <- "results/Seurat_output/cal07_NormalizedUMI_infected_2020-04-21.txt"
write.table(rbind(temp, temp2), file = out.file,
            row.names = FALSE, sep = "\t")



## Do tSNE of mock cells ----

so_mock <- subset(so, subset = Library == "Mock" & InfectedStatus == "NotInfected")


#Do SCTransform normalization
#NOTE: this will remove more genes because they aren't expressed in 
#just the mock cells

so_mock
so_mock <- SCTransform(so_mock, return.only.var.genes = FALSE)
so_mock

tail(rownames(so_mock))

#output normalized UMI values

temp <- as.matrix(so_mock[["SCT"]]@counts)
#wait to add seurat clusters


#remove viral genes redo PCA, and tSNE for mock cells ----
newNum <- nrow(so_mock)-7
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

out.file <- "results/Seurat_output/cal07_NormalizedUMI_Mock_2020-04-21.txt"

temp2 <- t(so_mock[[]][,c("seurat_clusters","NS_status")]) %>% as.data.frame() 
all.equal(colnames(temp), colnames(temp2))

rbind(as.data.frame(temp),temp2) %>% tibble::rownames_to_column(var = "gene") %>%
  write.table(file = out.file,
            row.names = FALSE, sep = "\t")



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


#Do tSNE plots for pres/abs each gene ----

for(i in c("PB2","PB1","PA","HA","NP","NA","M","NS")){
  plot <- DimPlot(so_infected, reduction = "tsne",group.by = paste0(i,"_status"),pt.size = 1)
  x11(width = 8, height = 7)
  plot + theme(title = element_text(size = rel(1.5)), plot.title = element_text(size = rel(2)), legend.text = element_text(size = rel(1.25)))
  out.file <- paste0("results/Seurat_output/cal07_Seurat_tSNE_Infected_",i,"_2020-04-20.jpeg")
  ggsave(out.file)
}



#Check a few other genes

plot <- FeaturePlot(so_infected,reduction = "tsne", features = "SLFN5",pt.size = 1)
x11(width = 8, height = 7)
plot + theme(title = element_text(size = rel(1.5)),  legend.text = element_text(size = rel(1.25)))
# out.file <- "results/Seurat_output/cal07_Seurat_tSNE_Infected_ViralPct_2020-04-13.jpeg"
# ggsave(out.file)

plot <- FeaturePlot(so_infected,reduction = "tsne", features = "IFNAR2",pt.size = 1)
x11(width = 8, height = 7)
plot + theme(title = element_text(size = rel(1.5)),  legend.text = element_text(size = rel(1.25)))
# out.file <- "results/Seurat_output/cal07_Seurat_tSNE_Infected_ViralPct_2020-04-13.jpeg"
# ggsave(out.file)



#Make RidgePlots for Fig 7 ----

#find the default colors to reverse them:
#From https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette

library(ggridges)

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

ggplotColours(n = 3)
# "#F8766D" "#00BA38" "#619CFF"

#NOTE: plots dated 2020-04-27 or 2020-04-28 had the colors in order "#00BA38","#F8766D","#619CFF"
#plots dated 2020-05-12 had the colors in order "#619CFF","#00BA38","#F8766D"


table(so$Library)
#  Mock Bystander  Infected 
#  4314      3367      4850

#subset to only cells with >0 for IFIT2 ----

DefaultAssay(so) <- "RNA"

so.sm <- subset(so, subset = IFIT2 > 0, slot = "counts")
#reorder to Infected Bystander, Mock
so.sm$Library <- factor(so.sm$Library, levels = c("Infected","Bystander","Mock"))
table(so.sm$Library)
#Infected Bystander      Mock 
#     195       187        67 

temp <- FetchData(so.sm, "IFIT2", slot = "counts")

table(temp >= 3, so.sm$Library)
#        Infected Bystander Mock
# FALSE      162       183   67
# TRUE        33         4    0

sum(temp[so.sm$Library == "Mock",1] >= 3) / sum(so$Library == "Mock") * 100
#0
sum(temp[so.sm$Library == "Bystander",1] >= 3) / sum(so$Library == "Bystander") * 100
#0.1188001
sum(temp[so.sm$Library == "Infected",1] >= 3) / sum(so$Library == "Infected") * 100
#0.6804124


x11(5,5)
RidgePlot(so.sm, features = "IFIT2", group.by = "Library", assay = "RNA",
          same.y.lims = TRUE, cols = c("#619CFF","#00BA38","#F8766D"),
          log = TRUE, slot = "counts") +  
  #geom_vline(xintercept = mock95, size = 1) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  xlab("raw counts") + ylab("") + NoLegend() +
  theme(plot.title = element_text(hjust = 0.3))

#ggsave("results/Seurat_output/cal07_RidgePlot_IFIT2_logAxis_2020-04-27.pdf")
#old colors

temp <- GetAssayData(object = so.sm, slot = "counts", assay = "RNA") %>% 
           as.matrix %>% log10()

so.sm <- SetAssayData(object = so.sm, assay = "RNA", slot = "data", 
                      new.data = temp)

x11(5,5)
RidgePlot(so.sm, features = "IFIT2", group.by = "Library", assay = "RNA",
          same.y.lims = TRUE, cols = c("#619CFF","#00BA38","#F8766D"),
          log = FALSE, slot = "data") +  
  geom_vline(xintercept = 0.38, size = 1) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  xlab("log10(raw UMI counts)") + ylab("") + NoLegend() +
  theme(plot.title = element_text(hjust = 0.3))

#ggsave("results/Seurat_output/cal07_RidgePlot_IFIT2_logCounts_2020-04-28.pdf")
#add in % of cells >= 3 counts by hand in Acrobat and saved as 
#cal07_RidgePlot_IFIT2_logCounts_%_2020-04-28.pdf
ggsave("results/Seurat_output/cal07_RidgePlot_IFIT2_logCounts_2020-05-12.pdf")




#Try plotting as percentage of total host counts

temp <- GetAssayData(object = so.sm, slot = "counts", assay = "RNA") %>% 
  as.matrix()

#change to % of total counts per sample

temp2 <- sweep(temp, MARGIN = 2, STATS = so.sm$host_counts, FUN = "/")
temp2 <- log10(temp2 * 100)

so.sm <- SetAssayData(object = so.sm, assay = "RNA", slot = "scale.data", 
                      new.data = temp2)

#What 95pctile of mock?

temp <- FetchData(so.sm, "IFIT2", slot = "scale.data")

mock95 <- quantile(temp[so.sm$Library == "Mock",1], 0.95)

#find percentage of all cells over mock95

sum(temp[so.sm$Library == "Mock",1] >= mock95) / sum(so$Library == "Mock") * 100
#0.09272137
sum(temp[so.sm$Library == "Bystander",1] >= mock95) / sum(so$Library == "Bystander") * 100
#0.4752005
sum(temp[so.sm$Library == "Infected",1] >= mock95) / sum(so$Library == "Infected") * 100
#1.381443

#find percentage of non-zero cells over mock95

over95 <- data.frame(Library = so.sm$Library, gene = temp[,1]) %>% group_by(Library) %>%
  summarize(pct = mean(gene > mock95) * 100) %>% select(pct) %>% round(1) %>% as.data.frame()
data.frame(Species = levels(so.sm$Library), over95t = paste0(over95[,1], "%") )
#    Species over95t
# 1  Infected     34.4%
# 2 Bystander      8.6%
# 3      Mock      6%

x11(5,5)
RidgePlot(so.sm, features = "IFIT2", group.by = "Library", assay = "RNA",
          same.y.lims = TRUE, cols = c("#619CFF","#00BA38","#F8766D"),
          log = FALSE, slot = "scale.data") +  
  geom_vline(xintercept = mock95, size = 1) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  xlab("log10(percent of total host counts)") + ylab("") + NoLegend() +
  theme(plot.title = element_text(hjust = 0.3))

#ggsave("results/Seurat_output/cal07_RidgePlot_IFIT2_PctTotalHost_2020-04-28.pdf")
#add in % of cells >= 3 counts by hand in Acrobat and saved as 
#cal07_RidgePlot_IFIT2_PctTotalHost_%_2020-04-28.pdf




#Do the same for ISG15

DefaultAssay(so) <- "RNA"

so.sm <- subset(so, subset = ISG15 > 0, slot = "counts")
#reorder to Infected Bystander, Mock
so.sm$Library <- factor(so.sm$Library, levels = c("Infected","Bystander","Mock"))
table(so.sm$Library)
#Infected Bystander      Mock 
#    1158       766      1805 

temp <- FetchData(so.sm, "ISG15", slot = "counts")

table(temp >= 3, so.sm$Library)
#        Infected Bystander Mock
# FALSE     1060       714 1578
# TRUE        98        52  227

sum(temp[so.sm$Library == "Mock",1] >= 3) / sum(so$Library == "Mock") * 100
#5.261938
sum(temp[so.sm$Library == "Bystander",1] >= 3) / sum(so$Library == "Bystander") * 100
#1.544402
sum(temp[so.sm$Library == "Infected",1] >= 3) / sum(so$Library == "Infected") * 100
#2.020619


x11(5,5)
RidgePlot(so.sm, features = "ISG15", group.by = "Library", assay = "RNA",
          same.y.lims = TRUE, cols = c("#619CFF","#00BA38","#F8766D"),
          log = TRUE, slot = "counts") +  
  #geom_vline(xintercept = mock95, size = 1) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  xlab("raw counts") + ylab("") + NoLegend() +
  theme(plot.title = element_text(hjust = 0.3))

#ggsave("results/Seurat_output/cal07_RidgePlot_ISG15_logAxis_2020-04-27.pdf")


temp <- GetAssayData(object = so.sm, slot = "counts", assay = "RNA") %>% 
  as.matrix %>% log10()

so.sm <- SetAssayData(object = so.sm, assay = "RNA", slot = "data", 
                      new.data = temp)

x11(5,5)
RidgePlot(so.sm, features = "ISG15", group.by = "Library", assay = "RNA",
          same.y.lims = TRUE, cols = c("#619CFF","#00BA38","#F8766D"),
          log = FALSE, slot = "data") +  
  geom_vline(xintercept = 0.38, size = 1) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  xlab("log10(raw UMI counts)") + ylab("") + NoLegend() +
  theme(plot.title = element_text(hjust = 0.3))

#ggsave("results/Seurat_output/cal07_RidgePlot_ISG15_logCounts_2020-04-28.pdf")
#add in % of cells >= 3 counts by hand in Acrobat and saved as 
#cal07_RidgePlot_ISG15_logCounts_%_2020-04-28.pdf
ggsave("results/Seurat_output/cal07_RidgePlot_ISG15_logCounts_2020-05-12.pdf")

#Try plotting as percentage of total host counts

temp <- GetAssayData(object = so.sm, slot = "counts", assay = "RNA") %>% 
  as.matrix()

#change to % of total counts per sample

temp2 <- sweep(temp, MARGIN = 2, STATS = so.sm$host_counts, FUN = "/")
temp2 <- log10(temp2 * 100)

so.sm <- SetAssayData(object = so.sm, assay = "RNA", slot = "scale.data", 
                      new.data = temp2)

#What 95pctile of mock?

temp <- FetchData(so.sm, "ISG15", slot = "scale.data")

mock95 <- quantile(temp[so.sm$Library == "Mock",1], 0.95)

#find percentage of all cells over mock95

sum(temp[so.sm$Library == "Mock",1] >= mock95) / sum(so$Library == "Mock") * 100
#2.109411
sum(temp[so.sm$Library == "Bystander",1] >= mock95) / sum(so$Library == "Bystander") * 100
#0.7722008
sum(temp[so.sm$Library == "Infected",1] >= mock95) / sum(so$Library == "Infected") * 100
#2.350515

#find percentage of non-zero cells over mock95

over95 <- data.frame(Library = so.sm$Library, gene = temp[,1]) %>% group_by(Library) %>%
  summarize(pct = mean(gene > mock95) * 100) %>% select(pct) %>% round(1) %>% as.data.frame()
data.frame(Species = levels(so.sm$Library), over95t = paste0(over95[,1], "%") )
#    Species over95t
# 1  Infected     9.8%
# 2 Bystander      3.4%
# 3      Mock      5%

x11(5,5)
RidgePlot(so.sm, features = "ISG15", group.by = "Library", assay = "RNA",
          same.y.lims = TRUE, cols = c("#619CFF","#00BA38","#F8766D"),
          log = FALSE, slot = "scale.data") +  
  geom_vline(xintercept = mock95, size = 1) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  xlab("log10(percent of total host counts)") + ylab("") + NoLegend() +
  theme(plot.title = element_text(hjust = 0.3))

#ggsave("results/Seurat_output/cal07_RidgePlot_ISG15_PctTotalHost_2020-04-28.pdf")
#add in % of cells >= 3 counts by hand in Acrobat and saved as 
#cal07_RidgePlot_ISG15_PctTotalHost_%_2020-04-28.pdf


#Do the same for ZC3HAV1

DefaultAssay(so) <- "RNA"

so.sm <- subset(so, subset = ZC3HAV1 > 0, slot = "counts")
#reorder to Infected Bystander, Mock
so.sm$Library <- factor(so.sm$Library, levels = c("Infected","Bystander","Mock"))
table(so.sm$Library)
#Infected Bystander      Mock 
#     745       496       538 

temp <- FetchData(so.sm, "ZC3HAV1", slot = "counts")

table(temp >= 3, so.sm$Library)
#        Infected Bystander Mock
# FALSE      729       490  535
# TRUE        16         6    3

sum(temp[so.sm$Library == "Mock",1] >= 3) / sum(so$Library == "Mock") * 100
#0.06954103
sum(temp[so.sm$Library == "Bystander",1] >= 3) / sum(so$Library == "Bystander") * 100
#0.1782002
sum(temp[so.sm$Library == "Infected",1] >= 3) / sum(so$Library == "Infected") * 100
#0.3298969


x11(5,5)
RidgePlot(so.sm, features = "ZC3HAV1", group.by = "Library", assay = "RNA",
          same.y.lims = TRUE, cols = c("#619CFF","#00BA38","#F8766D"),
          log = TRUE, slot = "counts") +  
  #geom_vline(xintercept = mock95, size = 1) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  xlab("raw counts") + ylab("") + NoLegend() +
  theme(plot.title = element_text(hjust = 0.3))

#ggsave("results/Seurat_output/cal07_RidgePlot_ZC3HAV1_logAxis_2020-04-27.pdf")


temp <- GetAssayData(object = so.sm, slot = "counts", assay = "RNA") %>% 
  as.matrix %>% log10()

so.sm <- SetAssayData(object = so.sm, assay = "RNA", slot = "data", 
                      new.data = temp)

x11(5,5)
RidgePlot(so.sm, features = "ZC3HAV1", group.by = "Library", assay = "RNA",
          same.y.lims = TRUE, cols = c("#619CFF","#00BA38","#F8766D"),
          log = FALSE, slot = "data") +  
  geom_vline(xintercept = 0.38, size = 1) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  xlab("log10(raw UMI counts)") + ylab("") + NoLegend() +
  theme(plot.title = element_text(hjust = 0.3))

#ggsave("results/Seurat_output/cal07_RidgePlot_ZC3HAV1_logCounts_2020-04-28.pdf")
#add in % of cells >= 3 counts by hand in Acrobat and saved as 
#cal07_RidgePlot_ZC3HAV1_logCounts_%_2020-04-28.pdf


#Need to overide the automatic bandwidth selection which is too small
#This needs to be done in the call to geom_density_ridges(), which is
#a few layers deep within RidgePlot() and there isn't a way to pass it
#through. Instead, call debug on an internal function Seurat:::SingleExIPlot().

debugonce(Seurat:::SingleExIPlot)

#Run the RidgePlot() code below which will pop up the debugger for Seurat:::SingleExIPlot.
#Step through the code by hitting enter until you get to this line of code:

#geom <- list(geom_density_ridges(scale = 4), theme_ridges(), 
#scale_y_discrete(expand = c(0.01, 0)), scale_x_continuous(expand = c(0,0)))

# Hit enter to run it once, then manually enter this code:

#geom <- list(geom_density_ridges(scale = 4, bandwidth = 0.01), theme_ridges(), 
#scale_y_discrete(expand = c(0.01, 0)), scale_x_continuous(expand = c(0,0)))

#Then you can either step through the rest of the code by hitting enter or click on the
#"Execute the remainer of the current function" 


x11(5,5)
RidgePlot(so.sm, features = "ZC3HAV1", group.by = "Library", assay = "RNA",
          same.y.lims = TRUE, cols = c("#619CFF","#00BA38","#F8766D"),
          log = FALSE, slot = "data") +  
  geom_vline(xintercept = 0.38, size = 1) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  xlab("log10(raw UMI counts)") + ylab("") + NoLegend() +
  theme(plot.title = element_text(hjust = 0.3))

ggsave("results/Seurat_output/cal07_RidgePlot_ZC3HAV1_logCounts_2020-05-12.pdf")





#Try plotting as percentage of total host counts

temp <- GetAssayData(object = so.sm, slot = "counts", assay = "RNA") %>% 
  as.matrix()

#change to % of total counts per sample

temp2 <- sweep(temp, MARGIN = 2, STATS = so.sm$host_counts, FUN = "/")
temp2 <- log10(temp2 * 100)

so.sm <- SetAssayData(object = so.sm, assay = "RNA", slot = "scale.data", 
                      new.data = temp2)

#What 95pctile of mock?

temp <- FetchData(so.sm, "ZC3HAV1", slot = "scale.data")

mock95 <- quantile(temp[so.sm$Library == "Mock",1], 0.95)

#find percentage of all cells over mock95

sum(temp[so.sm$Library == "Mock",1] >= mock95) / sum(so$Library == "Mock") * 100
#0.6258693
sum(temp[so.sm$Library == "Bystander",1] >= mock95) / sum(so$Library == "Bystander") * 100
#0.1782002
sum(temp[so.sm$Library == "Infected",1] >= mock95) / sum(so$Library == "Infected") * 100
#1.896907

#find percentage of non-zero cells over mock95

over95 <- data.frame(Library = so.sm$Library, gene = temp[,1]) %>% group_by(Library) %>%
  summarize(pct = mean(gene > mock95) * 100) %>% select(pct) %>% round(1) %>% as.data.frame()
data.frame(Species = levels(so.sm$Library), over95t = paste0(over95[,1], "%") )
#    Species over95t
# 1  Infected     12.3%
# 2 Bystander      1.2%
# 3      Mock      5%

x11(5,5)
RidgePlot(so.sm, features = "ZC3HAV1", group.by = "Library", assay = "RNA",
          same.y.lims = TRUE, cols = c("#619CFF","#00BA38","#F8766D"),
          log = FALSE, slot = "scale.data") +  
  geom_vline(xintercept = mock95, size = 1) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  xlab("log10(percent of total host counts)") + ylab("") + NoLegend() +
  theme(plot.title = element_text(hjust = 0.3))

#ggsave("results/Seurat_output/cal07_RidgePlot_ZC3HAV1_PctTotalHost_2020-04-28.pdf")
#add in % of cells >= 3 counts by hand in Acrobat and saved as 
#cal07_RidgePlot_ZC3HAV1_PctTotalHost_%_2020-04-28.pdf


# Make tSNE plots for selected genes in Fig 7 ----

grep("^IFN", rownames(so), value = TRUE)
#"IFNLR1" "IFNGR1" "IFNL2"  "IFNL1"  "IFNAR2" "IFNAR1" "IFNGR2"

DefaultAssay(so_infected) <-"RNA"
temp <- FetchData(so_infected, "IFNL1", slot = "counts")[,1]

so_infected$IFNL1_pres <- ifelse(temp > 0, "Positive","Negative")

plot <- DimPlot(so_infected,reduction = "tsne", group.by = "IFNL1_pres",
                    pt.size = 1, cols = c("lightgrey","red"))
x11(width = 8, height = 7)
plot + theme(title = element_text(size = rel(1.5)),  legend.text = element_text(size = rel(1.25))) + labs(title = "")
out.file <- "results/Seurat_output/cal07_Seurat_tSNE_IFNL1_2020-04-28.jpeg"
ggsave(out.file)



#Make violin plots ----

temp <- FetchData(so_infected, "NEAT1", slot = "counts")

write.csv(table(temp > 0, so_infected$seurat_clusters), 
          file = "results/Seurat_output/cal07_VlnPlot_NEAT1_2020-05-13.csv",
          row.names = c("zero","non-zero"))

x11(width = 8, height = 5)
VlnPlot(so_infected, features = "NEAT1", assay = "SCT",
        slot = "counts", log = TRUE) + NoLegend() + xlab("cluster") + ylab("normalized count") +
  theme(text = element_text(size = 20), axis.text.x= element_text(size = 15),
        axis.text.y= element_text(size = 15))

ggsave("results/Seurat_output/cal07_VlnPlot_NEAT1_2020-05-13.jpeg")



temp <- FetchData(so_infected, "TRIM41", slot = "counts")

write.csv(table(temp > 0, so_infected$seurat_clusters), 
          file = "results/Seurat_output/cal07_VlnPlot_TRIM41_2020-05-13.csv",
          row.names = c("zero","non-zero"))

x11(width = 8, height = 5)
VlnPlot(so_infected, features = "TRIM41", assay = "SCT",
        slot = "counts", log = TRUE) + NoLegend() + xlab("cluster") + ylab("normalized count") +
  theme(text = element_text(size = 20), axis.text.x= element_text(size = 15),
        axis.text.y= element_text(size = 15))

ggsave("results/Seurat_output/cal07_VlnPlot_TRIM41_2020-05-13.jpeg")



temp <- FetchData(so_infected, "TRIM28", slot = "counts")

write.csv(table(temp > 0, so_infected$seurat_clusters), 
          file = "results/Seurat_output/cal07_VlnPlot_TRIM28_2020-05-13.csv",
          row.names = c("zero","non-zero"))

x11(width = 8, height = 5)
VlnPlot(so_infected, features = "TRIM28", assay = "SCT",
        slot = "counts", log = TRUE) + NoLegend() + xlab("cluster") + ylab("normalized count") +
  theme(text = element_text(size = 20), axis.text.x= element_text(size = 15),
        axis.text.y= element_text(size = 15))

ggsave("results/Seurat_output/cal07_VlnPlot_TRIM28_2020-05-13.jpeg")


save.image("results/Seurat_output/cal07.RData")

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
#   [1] magrittr_1.5                dplyr_0.8.5                
# [3] svglite_1.2.3               ggplot2_3.3.0              
# [5] Seurat_3.1.4                DropletUtils_1.6.1         
# [7] SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.1
# [9] DelayedArray_0.12.2         BiocParallel_1.20.1        
# [11] matrixStats_0.56.0          Biobase_2.46.0             
# [13] GenomicRanges_1.38.0        GenomeInfoDb_1.22.1        
# [15] IRanges_2.20.2              S4Vectors_0.24.3           
# [17] BiocGenerics_0.32.0         simpleSingleCell_1.10.1    
# 
# loaded via a namespace (and not attached):
#   [1] systemfonts_0.1.1      sn_1.6-1               plyr_1.8.6            
# [4] igraph_1.2.5           lazyeval_0.2.2         splines_3.6.3         
# [7] listenv_0.8.0          TH.data_1.0-10         digest_0.6.25         
# [10] htmltools_0.4.0        gdata_2.18.0           fansi_0.4.1           
# [13] cluster_2.1.0          ROCR_1.0-7             limma_3.42.2          
# [16] globals_0.12.5         R.utils_2.9.2          sandwich_2.5-1        
# [19] colorspace_1.4-1       rappdirs_0.3.1         ggrepel_0.8.2         
# [22] xfun_0.12              crayon_1.3.4           RCurl_1.98-1.1        
# [25] jsonlite_1.6.1         survival_3.1-8         zoo_1.8-7             
# [28] ape_5.3                glue_1.4.0             gtable_0.3.0          
# [31] zlibbioc_1.32.0        XVector_0.26.0         leiden_0.3.3          
# [34] Rhdf5lib_1.8.0         future.apply_1.4.0     HDF5Array_1.14.3      
# [37] scales_1.1.0           mvtnorm_1.1-0          edgeR_3.28.1          
# [40] bibtex_0.4.2.2         Rcpp_1.0.4.6           metap_1.3             
# [43] plotrix_3.7-7          viridisLite_0.3.0      reticulate_1.15       
# [46] dqrng_0.2.1            rsvd_1.0.3             tsne_0.1-3            
# [49] htmlwidgets_1.5.1      httr_1.4.1             gplots_3.0.3          
# [52] RColorBrewer_1.1-2     TFisher_0.2.0          ellipsis_0.3.0        
# [55] ica_1.0-2              pkgconfig_2.0.3        R.methodsS3_1.8.0     
# [58] farver_2.0.3           uwot_0.1.8             locfit_1.5-9.4        
# [61] tidyselect_1.0.0       labeling_0.3           rlang_0.4.5           
# [64] reshape2_1.4.3         munsell_0.5.0          tools_3.6.3           
# [67] cli_2.0.2              ggridges_0.5.2         stringr_1.4.0         
# [70] npsurv_0.4-0           fitdistrplus_1.0-14    caTools_1.18.0        
# [73] purrr_0.3.3            RANN_2.6.1             pbapply_1.4-2         
# [76] future_1.16.0          nlme_3.1-145           R.oo_1.23.0           
# [79] compiler_3.6.3         rstudioapi_0.11        plotly_4.9.2.1        
# [82] png_0.1-7              lsei_1.2-0             tibble_3.0.0          
# [85] stringi_1.4.6          gdtools_0.2.2          lattice_0.20-38       
# [88] Matrix_1.2-18          multtest_2.42.0        vctrs_0.2.4           
# [91] mutoss_0.1-12          pillar_1.4.3           lifecycle_0.2.0       
# [94] Rdpack_0.11-1          lmtest_0.9-37          RcppAnnoy_0.0.16      
# [97] data.table_1.12.8      cowplot_1.0.0          bitops_1.0-6          
# [100] irlba_2.3.3            gbRd_0.4-11            patchwork_1.0.0       
# [103] R6_2.4.1               KernSmooth_2.23-16     gridExtra_2.3         
# [106] codetools_0.2-16       MASS_7.3-51.5          gtools_3.8.2          
# [109] assertthat_0.2.1       rhdf5_2.30.1           withr_2.1.2           
# [112] sctransform_0.2.1      mnormt_1.5-6           multcomp_1.4-12       
# [115] GenomeInfoDbData_1.2.2 grid_3.6.3             tidyr_1.0.2           
# [118] Rtsne_0.15             numDeriv_2016.8-1.1    tinytex_0.21     