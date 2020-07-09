#! /usr/bin/env Rscript

##get command line
args <- commandArgs(TRUE)
out <- args[4]

##load libraries
library(singleCellTK)
library(ccRemover)

##load object
SCE <- createSCE(assayFile = args[1],annotFile = args[2],featureFile = args[3])

##extract normailed/log2 test dataframe
logcounts_temp <- logcounts(SCE)

##center data by genes
mean_gene_exp <- rowMeans(logcounts_temp)
data_cen <- logcounts_temp - mean_gene_exp
remove(logcounts_temp)

##generate cell cycle gene index for human ensembl
gene_names <- rownames(data_cen)
cell_cycle_gene_indices <- gene_indexer(gene_names, species = "human",name_type = "ensembl" )

##generate final true/false vector
if_cc <- rep(FALSE,nrow(data_cen))
if_cc[cell_cycle_gene_indices] <- TRUE

##join results into list
dat <- list(x=data_cen, if_cc=if_cc)
remove(data_cen)

##run ccRemover
xhat <- ccRemover(dat, bar=FALSE)

##adjust matrix
xhat <- xhat + mean_gene_exp

##output new dataframe
write.table(xhat,file=out,sep="\t")

q()
