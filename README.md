# SingleCell
### 10X Chromium single cell pipeline scripts
#### For use in analysis of Influenza infected human lung cells

## Analysis Overview:
These scripts use the SimpleSingleCell and Seurat packages to perform various analyses downstream of Cell Ranger on a 10X Chromium Single Cell experiment (with Influenza infected human lung cells).  Steps performed include quality filtering, normalization, annotation, dimensional reduction (tSNE), and differential gene expression.

### requirements:
The scripts require either R or Perl to be installed and made available from the command line (with versions used during testing included in parentheses):
* R (v3.6.0)
* Perl (v5.24.2)

The following R packages are also required (with versions used during testing included in parentheses):
* simpleSingleCell (v1.8.0)
* EnsDb.Hsapiens.v86 (v2.99.0)
* DropletUtils (v1.4.1)
* scater (v1.12.2)
* scran (v1.12.1)
* BiocSingular (v1.0.0)
* Seurat (v3.0.2)
* ggplot2 (v3.1.1)
* svglite (v1.2.2)
* dplyr (v0.8.1)
* magrittr (v1.5)
* NBID (v0.1.1)

### preliminary steps: 
Raw reads should be demultiplexed and mapped to a host/virus hybrid reference using the 10X Chromium Cell Ranger software package:
https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

#### Hybrid Reference: 
-1- Download and unarchive the appropriate 10X single cell human reference (e.g. GRCh38-3.0.0):
https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

-2- Combine the virus genome sequence fasta file with the human genome sequence fasta file:
> e.g. cat GRCh38-3.0.0.FASTA CAL09.FASTA >GRCh38_CAL09.FASTA

-3- Add virus annotation to the human reference GTF annotation file.  This may require manual editing of the added virus rows to conform to GTF format standards (e.g. each virus gene should have gene, transcript, and exon rows; all rows should have unique ID tags).

-4- run the Cell Ranger Count function separately on each experimental library.

-5- run the Cell Ranger Aggr function on all experimental libraries to combine them into a single matrix.

### first script: Run_simpleSingleCell_GenerateObject.R
generate the filtered/normalized/annotated matrix from the raw Cell Ranger matrix output. This script should be run twice, the first time to generate a filtered/normalized matrix in text (tab-delimited) format and the unfiltered CellID list.  These tables are used to generate the metadata/factors that can then be imported by running the script a second time to produce the final, filtered/normalized/annotated matrix.

#### running the first script:
run the script from the command line using the following options:
input1: path to raw Cell Ranger matrix files
input2: output base file name

> e.g. Rscript Run_simpleSingleCell_GenerateObject.R Path/to/Raw/CellRanger/Matrix OutputBaseName

The script performs preliminary quality filtering and count normalization using the SimpleSingleCell package (primarily Scran functions), resulting in the following outputs:
output1: filtered, normalized matrix file (tab-delimited text file)
output2: raw matrix Cell ID list

#### generating the metadata table:
generate the required cell annotations using the 
