#! /usr/bin/env Rscript

##load first package
library(simpleSingleCell)

##get command line
args <- commandArgs(TRUE)
sce <- readRDS(args[1])   ###sce: input object file. reads the single cell experiment object
procs <- args[2]          ###procs: number of processors to use
mode <- args[3]           ###mode: 1=ANOVA-like; 2=pairwise (pick two groups)
fact <- args[4]           ###fact: factor to test
out <- args[5]            ###out: output base name (multiple files are output)
groups <- args[6]         ###groups: two group names to use for mode 2 (pairwise) analysis. comma-delimited, leave empty for two factor analysis

###defaults/settings
combine_background = 0   ###for pairwise individual virus gene comparisons: combine "absent" and "background" groups into "false" and checks for "infected" status
if (isEmpty(groups)){
  combine_background = 1
  groups <- "NA"
  print(groups)
  print("here!")
}else{
  groups <- unlist(strsplit(groups, ","))
}

##factor list
factors <- c("library", "TotalVirus", "TotalPB2", "TotalPB1", "TotalPA", "TotalHA", "TotalNP", "TotalNA", "TotalM2", "TotalM1", "TotalNS1", "TotalNEP", "StatusInfected", "StatusPB2", "StatusPB1", "StatusPA", "StatusHA", "StatusNP", "StatusNA", "StatusM2", "StatusM1", "StatusNS1", "StatusNEP", "MissingVirusGenes", "NumMissingVirusGenes", "Cluster")

##get scran size factors
scranFactors <- sizeFactors(sce)
##get cell metadata to use as factors
colFactors <- colData(sce)
##get count matrix
counts <- as.matrix(assay(sce,"counts"))

##main package
library(NBID)

##mode 1 or mode 2 for each factor
if (mode == 1){
  ##run diff expression
  result = DEUsingNBID(counts,colFactors[[fact]],ncore = procs,sizeFactor = scranFactors)
}else{
  if (fact == "library"){
    ##get index and subset by groups
    index = (colFactors$library == groups[1] | colFactors$library == groups[2])
    counts_sub = counts[, index]
    library_sub = colFactors$library[index]
    library_sub <- factor(library_sub)
    scranFactors_sub = scranFactors[index]
    result = DEUsingNBID(counts_sub,library_sub,ncore = procs,sizeFactor = scranFactors_sub)
  }
  ##MissingVirusGenes factor was calculated only for infected cells, so all other cells have a meaningless zero
  if (fact == "MissingVirusGenes"){
    ##special double test for both infected and the missing gene binary status (true/false)
    if (combine_background == 1){
      index = (colFactors$StatusInfected == "infected")
      counts_sub = counts[, index]
      MissingVirusGenes_sub = colFactors$MissingVirusGenes[index]
      scranFactors_sub = scranFactors[index]
      index_sub = (MissingVirusGenes_sub == 1)
      result = DEUsingNBID(counts_sub,index_sub,ncore = procs,sizeFactor = scranFactors_sub)
    }else{
      ##cannot perform this version of the test due to fake zeroes
      print("This test version has meaningless zero values!")
      q()
    }
  }
  if (fact == "StatusNS1"){
    ##special double test for both infected and the gene present binary status (true/false)
    if (combine_background == 1){
      index = (colFactors$StatusInfected == "infected")
      counts_sub = counts[, index]
      statusNS1_sub = colFactors$StatusNS1[index]
      scranFactors_sub = scranFactors[index]
      index_sub = (statusNS1_sub == "Present")
      result = DEUsingNBID(counts_sub,index_sub,ncore = procs,sizeFactor = scranFactors_sub)
    }else{
      index = (colFactors$StatusNS1 == groups[1] | colFactors$StatusNS1 == groups[2])
      counts_sub = counts[, index]
      statusNS1_sub = colFactors$StatusNS1[index]
      statusNS1_sub <- factor(statusNS1_sub)
      scranFactors_sub = scranFactors[index]
      result = DEUsingNBID(counts_sub,statusNS1_sub,ncore = procs,sizeFactor = scranFactors_sub)
    }
  }
  if (fact == "StatusNA"){
    ##special double test for both infected and the gene present binary status (true/false)
    if (combine_background == 1){
      index = (colFactors$StatusInfected == "infected")
      counts_sub = counts[, index]
      status_sub = colFactors$StatusNA[index]
      scranFactors_sub = scranFactors[index]
      index_sub = (status_sub == "Present")
      result = DEUsingNBID(counts_sub,index_sub,ncore = procs,sizeFactor = scranFactors_sub)
    }else{
      index = (colFactors$StatusNA == groups[1] | colFactors$StatusNA == groups[2])
      counts_sub = counts[, index]
      statusNA_sub = colFactors$StatusNA[index]
      statusNA_sub <- factor(statusNA_sub)
      scranFactors_sub = scranFactors[index]
      result = DEUsingNBID(counts_sub,statusNA_sub,ncore = procs,sizeFactor = scranFactors_sub)
    }
  }
  if (fact == "StatusPB2"){
    ##special double test for both infected and the gene present binary status (true/false)
    if (combine_background == 1){
      index = (colFactors$StatusInfected == "Infected")
      counts_sub = counts[, index]
      status_sub = colFactors$StatusPB2[index]
      scranFactors_sub = scranFactors[index]
      index_sub = (status_sub == "Present")
      result = DEUsingNBID(counts_sub,index_sub,ncore = procs,sizeFactor = scranFactors_sub)
    }else{
      index = (colFactors$StatusPB2 == groups[1] | colFactors$StatusPB2 == groups[2])
      counts_sub = counts[, index]
      status_sub = colFactors$StatusPB2[index]
      status_sub <- factor(status_sub)
      scranFactors_sub = scranFactors[index]
      result = DEUsingNBID(counts_sub,status_sub,ncore = procs,sizeFactor = scranFactors_sub)
    }
  }
}

##generate adjusted p-values (FDR) and print table
FDR = p.adjust(result[, "pvalue"], method = "BH")
result <- cbind(FDR, result)
write.table(result,file = out,sep = "\t")

q()
