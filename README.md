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

-4- run the Cell Ranger mkgtf function to filter the the new GTF file.

-5- run the Cell Ranger mkref function to index and prepare the new hybrid virus/host reference for use with Cell Ranger.


### first script: Run_simpleSingleCell_GenerateObject.R
generate the filtered/normalized/annotated matrix from the raw Cell Ranger matrix output. This script should be run twice, the first time to generate a filtered/normalized matrix in text (tab-delimited) format and the unfiltered CellID list.  These tables are used to generate the metadata/factors that can then be imported by running the script a second time to produce the final, filtered/normalized/annotated matrix.

#### running the first script, first time:
run the script from the command line once using the following options:
- input1: path to raw Cell Ranger matrix files
- input2: output base file name

> e.g. Rscript Run_simpleSingleCell_GenerateObject.R Path/to/Raw/CellRanger/Matrix OutputBaseName

The script performs preliminary quality filtering and count normalization using the SimpleSingleCell package (primarily Scran functions), resulting in the following outputs:
- output1: filtered, normalized matrix file (tab-delimited text file)
- output2: raw matrix Cell ID list

#### generating the metadata table:
generate the required cell annotations using the filtered, normalized matrix and combine with the raw matrix CellIDs to import back into the matrix object for downstream analysis (see example metadata table: MetaDataTable_Example.txt).

Required Factors/Cell annotations: Library, CellCycle, StatusInfected, StatusPB2, StatusPB1, StatusPA, StatusHA, StatusNP, StatusNA, StatusM, StatusNS, NumVirusGenes, ClusterID, Doublets

Optional Cell Annotations: TotalVirus, TotalPB2, TotalPB1, TotalPA, TotalHA, TotalNP, TotalNA, TotalM, TotalNS, AnyMissingVirusGenes, AnySingleMissingVirusGene

- Library: experimental library (e.g. treatment).
- CellCycle: cell cycle stage, as determined by the Cyclone tool in the Scran package (see Run_simpleSingleCell_Scran_CellCycle.R script).
- StatusInfected: cell infected status based on virus molecular count threshold (e.g. virus count versus expected virus background).
- Status\[VirusGene\]: virus gene presence/absence based on virus gene molecular count threshold (see StatusInfected).
- NumVirusGenes: number of virus genes present (i.e. expressed) in cell.
- Doublets: binary tag (i.e. 1 or 0) for droplets that are considered doublets for filtering purposes.

#### running the first script, second time:
Run the first script a second time using the following options:
- input1: path to raw Cell Ranger matrix files
- input2: output base file name
- input3: metadata table file name

The script performs preliminary quality filtering, count normalization, and matrix annotation using the SimpleSingleCell package (primarily Scran functions), resulting in the following outputs:
- output1: filtered, normalized, and annotated matrix file (RDS format) suitable for downstream analysis

### Seurat script: Run_Seurat_tSNEs_MASS.R
script to run Seurat and MAST functions in order to generate tSNE plots and DGE tables.

#### input1: filtered, normalized, and annotated matrix file (RDS format)

#### output1: high quality tSNE plots for the following:
     - all cells: Library, Cell Cycle, Total Virus (proportion)
     - infected cells: Library, Cell Cycle, Total Virus (proportion)
    
#### output2: DGE lists
     - all cells: Library, Seurat clusters
     -infected cells: Seurat clusters, StatusPb2, StatusPb1, StatusPA, StatusHA, StatusNP, StatusNA, StatusM, StatusNS


### Perl Script to Combine DGE tables: CompareLists_DGE.pl
combine differential gene expression tables (i.e. DGE lists) from different tools (i.e. NBID and MAST) or from different factors (e.g. StatusPB2, StatusPB1, and StatusPA) using a minimum intersection value. This allows results to be easily compared or for more robust differentially expressed genes to be selected (i.e. those genes called by both NBID and MAST).

#### inputs:
     -d: path to a directory containing DGE tables to be combined/compared.  Files must have the .tsv or .txt extention and have the following file name format: FileBaseName_\[MAST/NBID/Combined\]_\[FactorName\].[tsv\txt]
     -o: output file name
     -x: minimum intersection value
      
 #### outputs:
     -a tab-delimited table containing the ID column and the FDR/p-value and log2FC columns for each input file. 

> e.g. perl CompareLists_DGE.pl -d DGE/input/files/ -o OutputFileName.tsv -x 2

