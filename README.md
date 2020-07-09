# SingleCell
### 10X Chromium single cell scripts
#### Published as [Single cell heterogeneity in influenza A virus gene expression shapes the innate antiviral response to infection](https://doi.org/10.1371/journal.ppat.1008671), Sun et al. 2020 PLOS Pathogens

## Analysis Overview:
These scripts use the SimpleSingleCell and Seurat packages to perform various analyses downstream of Cell Ranger on a 10X Chromium Single Cell experiment (with Influenza infected human lung cells).  Steps performed include quality filtering, normalization, annotation, dimensional reduction (tSNE), and differential gene expression.

### requirements:
The scripts require R to be installed and made available from the command line (with versions used during testing included in parentheses):
* R (v3.6.3)


The following R packages are also required (with versions used during testing included in parentheses):
* simpleSingleCell (v1.10.1)
* EnsDb.Hsapiens.v86 (v2.99.0)
* DropletUtils (v1.6.1)
* scater (v1.14.6)
* scran (v1.14.6)
* BiocSingular (v1.2.2)
* Seurat (v3.1.4)
* ggplot2 (v3.3.0)
* svglite (v1.2.3)
* dplyr (v0.8.5)
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

-4- run the Cell Ranger mkgtf function to filter the the new GTF file.

-5- run the Cell Ranger mkref function to index and prepare the new hybrid virus/host reference for use with Cell Ranger.


### first script: Run_filter_and_count_virus_*.R
Generate the filtered/annotated matrix from the raw Cell Ranger matrix output with viral percentages and infection status calls. These scripts were originally designed to be called from the command line but were modified to run interactively due to slight differences between the two cell lines, cal07 and perth09.

The script performs preliminary quality filtering (empty drops, cell cycle calling, filter cells by detected genes, filter genes by min cells, double calling) and then calculates viral read percentages to call infection status and viral gene presence/absence, resulting in the following outputs:
- output1: non-normalized SingleCellExperiment object to be read into the Seurat analysis script
- output2: metadata.tsv file for all cells including viral percentages and calls
- output3: normalized SingleCellExperiment object of only infected cells to be read into the NBID analysis script
- output4: variety of graphs for the publication


### Seurat script: Seurat_*_tSNEs_MAST.R
Script to run Seurat and MAST functions in order to generate tSNE plots and DGE tables.

#### input1: non-normalized SingleCellExperiment object (RDS format)

#### outputs
- output1: high quality tSNE plots for the following:
     - all cells: Library, Cell Cycle, Total Virus (proportion)
     - infected cells: Library, Cell Cycle, Total Virus (proportion)
     - mock cells: Library
 - output2: Normalized UMI count matrix for infected cells
 - output3: DGE lists from MAST test
 - output4: variety of graphs for the publication
   
    
### NBID script: NBID_NCore_VirGenes_*.sh
These shell scripts in turn call Run_simpleSingleCell_NBID_NCore_VirGenes.R to run the NBID differential expression testing using multiple cores on a computer cluster using a SLURM scheduler.

#### input1: Normalized SingleCellExperiment object of infected cells only (RDS format)

#### outputs: DGE lists comparing cells expressing or not expressing each viral gene

 
### Script to Combine DGE tables: CompareLists_DGE_*.R
combine differential gene expression tables (i.e. DGE lists) from different tools (i.e. NBID and MAST) or from different factors (e.g. StatusPB2, StatusPB1, and StatusPA) using a minimum intersection value. This allows results to be easily compared or for more robust differentially expressed genes to be selected (i.e. those genes called by both NBID and MAST).

#### inputs:
    - .txt files output from Seurat_*_tSNEs_MAST.R
    - .txt files output from NBID_NCore_VirGenes_*.sh
      
 #### outputs:
     - a tab-delimited table for each viral gene showing DE genes from both tests and whether each DE gene is unique to that viral gene.
     - a tab-delimited table with significant and unique DE gene counts for each viral gene.
     


