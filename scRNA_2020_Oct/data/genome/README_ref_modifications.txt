Used this procedure for salmon: https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

cd /home/groups/hpcbio/RNA-Seq/projects/ningwang/2020Sep-RNASeq/data/genome

#Get the references

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/223/135/GCF_000223135.1_CriGri_1.0/GCF_000223135.1_CriGri_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/223/135/GCF_000223135.1_CriGri_1.0/GCF_000223135.1_CriGri_1.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/223/135/GCF_000223135.1_CriGri_1.0/README_Cricetulus_griseus_annotation_release_104
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/223/135/GCF_000223135.1_CriGri_1.0/annotation_hashes.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/223/135/GCF_000223135.1_CriGri_1.0/GCF_000223135.1_CriGri_1.0_rna.fna.gz

gunzip *gz


#get the genome seqnames for the decoys list
grep "^>" GCF_000223135.1_CriGri_1.0_genomic.fna | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

#cat together the transcriptome and genome

cat GCF_000223135.1_CriGri_1.0_rna.fna GCF_000223135.1_CriGri_1.0_genomic.fna > gentrome.fa
gzip gentrome.fa


#modify gff file

srun --pty -p hpcbio bash
module load R/3.6.0-IGB-gcc-8.2.0
Rscript ../../src/gff2gtf_rtracklayer.R GCF_000223135.1_CriGri_1.0_genomic.gff GCF_000223135.1_CriGri_1.0_genomic__ModFromGFF.gtf putEGin_gene_id
exit


#The STAR reference was made with this script 

cd ../../src
sbatch star_index.sh


#The Salmon reference was made with this script 

sbatch salmon_index.sh
