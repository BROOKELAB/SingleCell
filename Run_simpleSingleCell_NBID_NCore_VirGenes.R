#! /usr/bin/env Rscript

##to install NBID
#install.packages("devtools")
#library(devtools)
#install_bitbucket("Wenan/NBID", ref="default", build_vignettes = T)

##packages installed
##simpleSingleCell_1.8.0
##DropletUtils_1.4.1
##NBID_0.1.1

##load packages
library(simpleSingleCell)
library(DropletUtils)
library(NBID)

##get command line
args <- commandArgs(TRUE)
sce.inf <- readRDS(args[1])   ###sce: input object file. reads the single cell experiment object
out <- args[2]
procs <- as.integer(args[3])          ###procs: number of processors to use 

##potential factors list: this list is not actually used below, it's here for conveniance only!
#factors2test <- c("Library", "CellCycle", "TotalVirus", "TotalPB2", "TotalPB1", "TotalPA", "TotalHA", "TotalNP", "TotalNA", "TotalM", "TotalNS", "StatusInfected", "StatusPB2", "StatusPB1", "StatusPA", "StatusHA", "StatusNP", "StatusNA", "StatusM", "StatusNS", "MissingVirusGenes", "SingleMissingVirusGene", "NumMissingVirusGenes", "ClusterID")

##get scran size factors
scranFactors <- sizeFactors(sce.inf)
##get cell metadata to use as factors
colFactors <- colData(sce.inf)
##get count matrix
counts <- as.matrix(assay(sce.inf,"counts"))

##individual missing virus genes in infected cells
##for infected cells on missing PB2
status_sub = colFactors$PB2_status
index_sub = (status_sub == "Present")
result = DEUsingNBID(counts,index_sub,sizeFactor = scranFactors,ncore = procs)
FDR = p.adjust(result[, "pvalue"], method = "BH")
result <- cbind(FDR, result)
out.file <- paste(out,"_Scran_NBID_InfectedCells_PB2.txt",sep = "")
write.table(result,file = out.file,sep = "\t")

##for infected cells on missing PB1
status_sub = colFactors$PB1_status
index_sub = (status_sub == "Present")
result = DEUsingNBID(counts,index_sub,sizeFactor = scranFactors,ncore = procs)
FDR = p.adjust(result[, "pvalue"], method = "BH")
result <- cbind(FDR, result)
out.file <- paste(out,"_Scran_NBID_InfectedCells_PB1.txt",sep = "")
write.table(result,file = out.file,sep = "\t")

##for infected cells on missing PA
status_sub = colFactors$PA_status
index_sub = (status_sub == "Present")
result = DEUsingNBID(counts,index_sub,sizeFactor = scranFactors,ncore = procs)
FDR = p.adjust(result[, "pvalue"], method = "BH")
result <- cbind(FDR, result)
out.file <- paste(out,"_Scran_NBID_InfectedCells_PA.txt",sep = "")
write.table(result,file = out.file,sep = "\t")

##for infected cells on missing HA
status_sub = colFactors$HA_status
index_sub = (status_sub == "Present")
result = DEUsingNBID(counts,index_sub,sizeFactor = scranFactors,ncore = procs)
FDR = p.adjust(result[, "pvalue"], method = "BH")
result <- cbind(FDR, result)
out.file <- paste(out,"_Scran_NBID_InfectedCells_HA.txt",sep = "")
write.table(result,file = out.file,sep = "\t")

##for infected cells on missing NP
status_sub = colFactors$NP_status
index_sub = (status_sub == "Present")
result = DEUsingNBID(counts,index_sub,sizeFactor = scranFactors,ncore = procs)
FDR = p.adjust(result[, "pvalue"], method = "BH")
result <- cbind(FDR, result)
out.file <- paste(out,"_Scran_NBID_InfectedCells_NP.txt",sep = "")
write.table(result,file = out.file,sep = "\t")

##for infected cells on missing NA
status_sub = colFactors$NA_status
index_sub = (status_sub == "Present")
result = DEUsingNBID(counts,index_sub,sizeFactor = scranFactors,ncore = procs)
FDR = p.adjust(result[, "pvalue"], method = "BH")
result <- cbind(FDR, result)
out.file <- paste(out,"_Scran_NBID_InfectedCells_NA.txt",sep = "")
write.table(result,file = out.file,sep = "\t")

##for infected cells on missing M
status_sub = colFactors$M_status
index_sub = (status_sub == "Present")
result = DEUsingNBID(counts,index_sub,sizeFactor = scranFactors,ncore = procs)
FDR = p.adjust(result[, "pvalue"], method = "BH")
result <- cbind(FDR, result)
out.file <- paste(out,"_Scran_NBID_InfectedCells_M.txt",sep = "")
write.table(result,file = out.file,sep = "\t")

##for infected cells on missing NS
status_sub = colFactors$NS_status
index_sub = (status_sub == "Present")
result = DEUsingNBID(counts,index_sub,sizeFactor = scranFactors,ncore = procs)
FDR = p.adjust(result[, "pvalue"], method = "BH")
result <- cbind(FDR, result)
out.file <- paste(out,"_Scran_NBID_InfectedCells_NS.txt",sep = "")
write.table(result,file = out.file,sep = "\t")



q()
