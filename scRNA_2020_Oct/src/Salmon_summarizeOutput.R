#Code to summarize Salmon output from quant.sf files and meta_info.json files

#run on worker node and in directory containing all Salmon output folders (one per sample):
# $ srun --pty bash
# $ module load R/3.5.0-IGB-gcc-4.9.4
# $ Rscript pathToFile/Salmon_summarizeOutput.R

#output is SalmonSummarizedOutput.RData file with R objects tx.all (contains abundances (TPM), counts and estimated lengths)
#and meta_info containing Salmon's output summary info including salmon version and parameters, read numbers input and mapped.
#See https://salmon.readthedocs.io/en/latest/file_formats.html#fileformats for information on Salmon's outputs

#If you also prefer the script to output tab-delimited text files as well, run the script like this:

# $ Rscript pathToFile/Salmon_summarizeOutput.R outputText


#R codes start here:

args <- commandArgs(TRUE)

library(tximport)
library(readr)
library(jsonlite)
library(ggplot2)

options(stringsAsFactors=F)


# Get the full paths to the quant.sf files in each subdirectory

files.quant <- dir(pattern = "quant.sf$", recursive = T, full.names = T)


# create wrapper around read_tsv function to automatically read in  the Name as character
# and everything as double format (needed if only integer values in first 1000 rows)

my_read_tsv <- function(x, ...) {
    read_tsv(x, col_types = cols(Name = "c", .default = "d"), ...)
}


tx.all <- tximport(files.quant, type = "salmon", txOut = TRUE, dropInfReps = TRUE, 
                   importer = my_read_tsv)
               
#add in cleaned sample names

if (length(grep("_L00", files.quant)) == length(files.quant)) {
    temp <- strsplit(files.quant, "_")
    quant.good <- sapply(temp, function(x) {
        Lpos <- grep("L00", x)
        paste(x[1:(Lpos-2)], collapse = "_")
    })
    } else {
    quant.good <- gsub("/quant.sf", "", files.quant)
    }

quant.good <- make.names(gsub("./", "", quant.good))

colnames(tx.all$abundance) <- quant.good
colnames(tx.all$counts) <- quant.good
colnames(tx.all$length) <- quant.good



#also read in and summarize the info in the meta_info.json: Get the full paths to the files in each subdirectory

files.info <- dir(pattern = "meta_info.json", recursive = T, full.names = T)


# read in the first file

temp <- fromJSON(files.info[1], flatten = T)

meta_info <- as.data.frame(temp[sapply(temp, length) == 1])

#Now do for the rest of the samples
#skip this step if there's only one sample
if(length(files.info)>1){
for (i in 2:length(files.info)) {
    meta_info <- rbind(meta_info, as.data.frame(fromJSON(files.info[i], flatten = T)[sapply(temp, length) == 1]))
}}

#add in cleaned sample names

if (length(grep("_L00", files.info)) == length(files.info)) {
    temp <- strsplit(files.info, "_")
    files.good <- sapply(temp, function(x) {
        Lpos <- grep("L00", x)
        paste(x[1:(Lpos-2)], collapse = "_")
    })
    } else {
    files.good <- gsub("/aux_info/meta_info.json", "", files.info)
    }

files.good <- make.names(gsub("./", "", files.good))

rownames(meta_info) <- files.good


#write out tx_all and meta_info objects

save(meta_info, tx.all, file = "SalmonSummarizedOutput.RData")

#also write out as tab-delimited text files if specified:

if ("outputText" %in% args) {
  write.table(data.frame(sample = rownames(meta_info), meta_info), file = "Salmon_MetaInfo.txt", 
              row.names = FALSE, sep = "\t")
  write.table(data.frame(transcript_id = rownames(tx.all$abundance), tx.all$abundance),
              file = "Salmon_Tx_TPMabundances.txt",row.names = FALSE, sep = "\t") 
  write.table(data.frame(transcript_id = rownames(tx.all$counts), tx.all$counts),
              file = "Salmon_Tx_NumReads.txt", row.names = FALSE, sep = "\t")
  write.table(data.frame(transcript_id = rownames(tx.all$length), tx.all$length),
              file = "Salmon_Tx_EffectiveLengths.txt", row.names = FALSE, sep = "\t")
}


#make pct mapped plot

meta_info$Sample <- factor(rownames(meta_info), levels = unique(rownames(meta_info)))

jpeg("SalmonPctMapped.jpeg", width = 10, height = 6, units = "in", res = 300, quality = 100, type = "cairo")
ggplot(data=meta_info, aes(x=Sample, y=percent_mapped)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  ggtitle("Percentage of reads mapped to transcriptome by Salmon")+ ylab("Percentage")+
  scale_y_continuous(breaks = seq(0,100,by = 10), limits = c(0,100)) 
dev.off()

