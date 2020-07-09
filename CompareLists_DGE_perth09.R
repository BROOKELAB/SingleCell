#### Compare MAST and NBID test results ####

library(dplyr)

fdr_min <- 0.01
out <- "results/CompareLists_DGE/2020-04-20_perth09_SigGenes_"

inMAST1 <- "results/Seurat_output/perth09_Seurat_MAST_DGElist_InfectedCells_"
inMAST2 <- "_2020-04-13.txt"

inNBID1 <- "results/NBID_output/2020-04-18-perth09_Scran_NBID_InfectedCells_"
inNBID2 <- ".txt"


allDGE <- data.frame(gene = c("PB2","PB1","PA","HA","NP","virNA","M","NS"), 
                     dgeMAST = NA, dgeNBID = NA, dgeBOTH = NA, dgeUNIQUE = NA)
rownames(allDGE) <- allDGE$gene
BOTHdge <- list(PB2=NA,PB1=NA,PA=NA,HA=NA,NP=NA, virNA = NA,M=NA,NS=NA)

for (i in allDGE$gene) {
  
  MASTdge <- read.delim(paste0(inMAST1, gsub("vir","",i), inMAST2))
  NBIDdge <- read.delim(paste0(inNBID1, gsub("vir","",i), inNBID2))
  #put gene/rownames in column; Seurat had replaced underscores with dashes,
  #so do the same for NBID results
  MASTdge<- tibble::rownames_to_column(MASTdge, var = "gene")
  NBIDdge$gene <- gsub("_","-", rownames(NBIDdge))
  
  if(sum(MASTdge$p_val_adj < fdr_min) == 0 | sum(NBIDdge$FDR < fdr_min) == 0){
    allDGE[i,"dgeBOTH"] <- 0
    allDGE[i,"dgeMAST"] <- sum(MASTdge$p_val_adj < fdr_min)
    allDGE[i,"dgeNBID"] <- sum(NBIDdge$FDR < fdr_min) 
  } else {
    MASTdge <- MASTdge[MASTdge$p_val_adj < fdr_min,]
    NBIDdge <- NBIDdge[NBIDdge$FDR < fdr_min,]
    allDGE[i,"dgeMAST"] <- nrow(MASTdge)
    allDGE[i,"dgeNBID"] <- nrow(NBIDdge)
    allDGE[i,"dgeBOTH"] <- length(intersect(MASTdge$gene, NBIDdge$gene))
    if (allDGE[i,"dgeBOTH"] > 0) {
      names(MASTdge) <- c("gene","pval_MAST","logFC_MAST", "pct.1","pct.2", "FDR_MAST" )
      names(NBIDdge) <- c("FDR_NBID","pval_NBID","LR","betaTRUE","dispFALSE","dispTRUE","logFC_NBID","gene")
      BOTHdge[[i]] <- inner_join(MASTdge, NBIDdge[,c(1,2,7,8)])
      }
  }
}

allDGE


# Now find how many genes unique to each viral gene

names(BOTHdge)
BOTHdge <- BOTHdge[allDGE$dgeBOTH > 0]


for (i in names(BOTHdge)){
  one <- BOTHdge[[i]]
  rest <- BOTHdge[-which(names(BOTHdge) == i)] %>% plyr::ldply(data.frame)
  one$unique <- !one$gene %in% rest$gene
  allDGE[i,"dgeUNIQUE"] <- sum(one$unique)
  if(sum(one$unique) > 0){
    one$others <- ""
    for (j in which(!one$unique)) {
      one$others[j] <- rest[rest$gene %in% one$gene[j],1] %>% 
        gsub("vir","",.) %>% paste(collapse = ",")
    }
  }
  write.table(one, file = paste0(out, gsub("vir","",i),".txt"),
              row.names = FALSE, sep = "\t")
}

allDGE

write.csv(allDGE, file = "results/CompareLists_DGE/2020-04-20_perth09_SigCounts.csv",
          row.names = FALSE)
