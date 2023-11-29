#!/bin/bash

# CHANGE TO DIRECTORY WITH BAM FILES AND RUN THIS SCRIPT FROM THERE
# RUN AS SUDO, IF STAR ALIGNMENT WAS ALSO RUN AS SUDO



########### DEFINE INPUTS (CHANGE AS NEEDED)

# Define the number of TEcount iterations to run in parallel
TEcount_processes=3

# Define gene/TE annotation files
gene_gtf="/Users/juanb/Documents/Bioinformatic_Tools/STAR_genome_indices/homo_sapiens_gencode_v33_plus_EBV/hs_gencode33_with_EBV.gtf"
TE_file="/Users/juanb/Documents/Bioinformatic_Tools/STAR_genome_indices/TElocal_rmsk_files/GRCh38_GENCODE_rmsk_TE.gtf.locInd"

# Define strandness of the sequencing library
TETranscripts_strandess='no'


########### GENE/TE QUANTIFICATION WITH TETRANSCRIPTS


# Make a directory to hold the counts tables.
mkdir counts

# Analyze files with TEcounts in batches of 'TEcount_processes.' Run in "multi" mode here.
ls *.bam | xargs -P $TEcount_processes -Ibam_file TElocal --BAM bam_file --GTF $gene_gtf --TE $TE_file --project bam_file --mode multi --sortByPos --stranded $TETranscripts_strandess

# Move all counts tables to their own folder.
mv *.cntTable counts
