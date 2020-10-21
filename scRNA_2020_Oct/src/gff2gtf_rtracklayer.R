# R script to convert NCBI's gff to gtf using rtracklayer package.
# Call from command line using:
# Rscript gff2gtf_rtracklayer.R existing.gff outputFileName.gtf putEGin_gene_id
# if specify putEGin_gene_id, the script will pull out the Entrez GeneIDs from
# the Dbxref column and put them into the attribute column under the standard "gene_id"
# Also, if no information in Dbxref column, use info in product column

# IsoformSwitchAnalyzeR expects all CDS and exons to have transcript_id's, so
# this script can also copy those over from the parent mRNA (txpt_id_from_parent).
# GFF3 format has the stop codon included in CDS whereas GTF does not, so the
# trim_stop_codon option adjusts the CDS ranges accordingly.  This is also important
# for IsoformSwitchAnalyzeR.

args <- commandArgs(TRUE)
if(any(!args[-(1:2)] %in% c("putEGin_gene_id", "txpt_id_from_parent", "trim_stop_codon"))){
  stop(paste("Unrecognized argument:", paste(args[-(1:2)][!args[-(1:2)] %in% c("putEGin_gene_id", "txpt_id_from_parent", "trim_stop_codon")])))
}

gff0 <- rtracklayer::import(args[1])

if("putEGin_gene_id" %in% args) {
  temp <- sapply(gff0$Dbxref, function(x) x[grep("GeneID", x)])
  #If didn't have GeneID, switch to NA
  temp2 <- sapply(temp, length)
  if(any(temp2 == 0))
    temp[temp2 == 0] <- NA
  #If have more than one GeneID, keep first one
  if(any(temp2 > 1))
    temp[temp2 > 1] <- lapply(temp[temp2 > 1], function(x) x[1])
  #Put GeneID in attribute named gene_id
  gff0$gene_id <- unlist(temp)
  #If no GeneID, put in value from product attribute                              
  if(sum(is.na(gff0$gene_id)) > 0)
    gff0$gene_id[is.na(gff0$gene_id)] <- gff0$product[is.na(gff0$gene_id)]
}

# setup for two subroutines below
if("txpt_id_from_parent" %in% args || "trim_stop_codon" %in% args){
  txpt_rows <- which(!is.na(gff0$transcript_id))
  transcripts <- gff0$transcript_id[txpt_rows]
  names(transcripts) <- gff0$ID[txpt_rows]
  
  # Change the parent column from list into vector
  parentcol <- rep(NA_character_, length(gff0))
  hasparent <- which(lengths(gff0$Parent) > 0)
  if(all(lengths(gff0$Parent[hasparent]) == 1)){
    parentcol[hasparent] <- unlist(gff0$Parent[hasparent])
  } else {
    # if there are multiple parents, pick the first that has a txpt
    parentcol[hasparent] <-
      sapply(gff0$Parent[hasparent],
             function(x){
               x <- x[x %in% names(transcripts)]
               if(length(x) == 0){
                 return(NA_character_)
               } else {
                 return(x[1])
               }
             } )
  }
}

# Propagate transcript_id to CDS and exons if necessary
if("txpt_id_from_parent" %in% args){
  cds_rows <- which(gff0$type == "CDS" & is.na(gff0$transcript_id)) # only CDS that need it filled in
  cds_rows <- cds_rows[which(parentcol[cds_rows] %in% names(transcripts))] # need parent with transcript_id
  if(length(cds_rows) > 0){
    gff0$transcript_id[cds_rows] <- transcripts[parentcol[cds_rows]]
  }
  
  exon_rows <- which(gff0$type == "exon" & is.na(gff0$transcript_id)) # exons that need to be filled in
  exon_rows <- exon_rows[which(parentcol[exon_rows] %in% names(transcripts))] # need parent with transcript_id
  if(length(exon_rows) > 0){
    gff0$transcript_id[exon_rows] <- transcripts[parentcol[exon_rows]]
  }
  
  # Fill in for those without transcript; gene_id if present, else product
  if(!is.null(gff0$gene_id)){
    tofill <- which(gff0$type %in% c("CDS", "exon") & is.na(gff0$transcript_id))
    gff0$transcript_id[tofill] <- gff0$gene_id[tofill]
  }
  if(!is.null(gff0$product)){
    tofill <- which(gff0$type %in% c("CDS", "exon") & is.na(gff0$transcript_id))
    gff0$transcript_id[tofill] <- gff0$product[tofill]
  }
}

# remove stop codon from CDS to comply with GTF format
if("trim_stop_codon" %in% args){
  # for exon-exon junctions in the middle of a stop codon
  removerows <- integer(0)
  
  # fill in for top strand
  CDS_top <- which(as.logical(gff0$type == "CDS" & BiocGenerics::strand(gff0) == "+"))
  CDS_top_split <- split(CDS_top, parentcol[CDS_top])
  modifyrows <- integer(length(CDS_top_split))
  tosubtract <- integer(length(CDS_top_split))
  for(i in seq_along(CDS_top_split)){
    cds <- CDS_top_split[[i]]
    last <- cds[which.max(BiocGenerics::start(gff0[cds]))]
    if(BiocGenerics::width(gff0[last]) > 3){
      # whole stop codon is in this exon (vast majority of cases)
      modifyrows[i] <- last
      tosubtract[i] <- 3L
    } else {
      # stop codon goes into the previous exon
      tosubtract[i] <- 3L - BiocGenerics::width(gff0[last]) # remaining nucleotides to remove
      removerows <- c(removerows, last)
      cds <- cds[cds != last]
      modifyrows[i] <- cds[which.max(start(gff0[cds]))]
    }
  }
  BiocGenerics::end(gff0[modifyrows]) <- BiocGenerics::end(gff0[modifyrows]) - tosubtract
  
  # fill in for bottom strand
  CDS_bot <- which(as.logical(gff0$type == "CDS" & BiocGenerics::strand(gff0) == "-"))
  CDS_bot_split <- split(CDS_bot, parentcol[CDS_bot])
  modifyrows <- integer(length(CDS_bot_split))
  tosubtract <- integer(length(CDS_bot_split))
  for(i in seq_along(CDS_bot_split)){
    cds <- CDS_bot_split[[i]]
    last <- cds[which.min(BiocGenerics::start(gff0[cds]))]
    if(BiocGenerics::width(gff0[last]) > 3){
      # whole stop codon is in this exon (vast majority of cases)
      modifyrows[i] <- last
      tosubtract[i] <- 3L
    } else {
      # stop codon goes into the previous exon
      tosubtract[i] <- 3L - BiocGenerics::width(gff0[last]) # remaining nucleotides to remove
      removerows <- c(removerows, last)
      cds <- cds[cds != last]
      modifyrows[i] <- cds[which.min(BiocGenerics::start(gff0[cds]))]
    }
  }
  BiocGenerics::start(gff0[modifyrows]) <- BiocGenerics::start(gff0[modifyrows]) + tosubtract
  
  if(length(removerows) > 0) gff0 <- gff0[-removerows] # CDS that only contain stop codon
}

rtracklayer::export(gff0, args[2], format = "gtf")
