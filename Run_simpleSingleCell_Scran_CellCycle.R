#! /usr/bin/env Rscript

##get command line
args <- commandArgs(TRUE)
indir <- args[1]
out <- args[2]


##load first libraries
library(simpleSingleCell)
library(DropletUtils)

##load object
sce <- read10xCounts(indir, col.names=TRUE)

##create cell cycle scores
library(scran)
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
cellcycles <- cyclone(sce, pairs=hs.pairs)

##output new dataframe
write.table(cellcycles,file=out,sep="\t")

q()
