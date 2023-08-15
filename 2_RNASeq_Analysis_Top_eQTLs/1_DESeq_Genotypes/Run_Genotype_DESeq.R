# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/1_DESeq_Genotypes/Run_Genotype_DESeq_functions.R") 


# Load libraries
library(DESeq2) # For differential expression 
library(BiocParallel) # Used by DESeq2 for parallelization
  register(MulticoreParam(6)) # Use six cores
library(BEDMatrix) # to load genotypes







# DESEQ REF VS ALT ALLELES
  
  
# Define output directories
dir.output.DESeq <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/1_DESeq_Genotypes/DESeq_Results/'  

  
  
  
  
  
  

# LOAD COUNTS AND KNOWN COVARIATES

# Load filtered counts
counts.YRI <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_eQTL/All_counts_YRI_86_filtered.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')
counts.EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_eQTL/All_counts_EUR_358_filtered.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Load the covariates data
covariates_table.YRI  <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_eQTL/COVARIATES_YRI.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
covariates_table.YRI <- covariates_table.YRI[colnames(counts.YRI), ]      

covariates_table.EUR  <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_eQTL/COVARIATES_EUR.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
covariates_table.EUR <- covariates_table.EUR[colnames(counts.EUR), ]   

# Load structural variant PCA covariates
SV.PCs <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/Population_Structure_Analysis/COVARIATE_Structural_Variant_Population_Structure_PCA.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # Assign unique names to PCs (since SNV genotype SVs are already labeled PC1-PC20)
    colnames(SV.PCs) <- paste('SV_PC', 1:ncol(SV.PCs), sep = '')

# Add SV PCs 1-3 to the covariates tables
covariates_table.YRI <- cbind(covariates_table.YRI, 
                              SV.PCs[colnames(counts.YRI), c('SV_PC1', 'SV_PC2', 'SV_PC3')]
                              )
covariates_table.EUR <- cbind(covariates_table.EUR, 
                              SV.PCs[colnames(counts.EUR), c('SV_PC1', 'SV_PC2', 'SV_PC3')]
                              )






# DEFINE GENOTYPES FOR TOP SNPS

# Specify the genotype BED files
plink_file_path.YRI <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_YRI_Genotypes/Filtered_YRI_BOTH_SEXES_chr1-22.recode.plink.bed'
plink_file_path.EUR <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink.bed'

# Make R object to stream snp genotypes into memory (rows are samples and columns are snps). # NOTE: If used, MT snps will need to be assigned a snpid, since they dont have rsid.
binary_genotypes.YRI <- BEDMatrix(path = plink_file_path.YRI, simple_names = T)
binary_genotypes.EUR <- BEDMatrix(path = plink_file_path.EUR, simple_names = T)

# Define the snps of interest
top.snps.YRI <- c('rs2176598', 'rs9271379') 
top.snps.EUR <- c('rs11635336', 'rs9271894', 'rs1061810', 'rs112581165', 'rs72691418', 'rs9270493') # NOTE: rs9270493 is for RNF5

# Extract genotypes for snps to be tested
top.genotypes.YRI <- as.data.frame(binary_genotypes.YRI[colnames(counts.YRI), top.snps.YRI])
top.genotypes.EUR <- as.data.frame(binary_genotypes.EUR[colnames(counts.EUR), top.snps.EUR])

# Add SNPs to covariates table
covariates_table.YRI <- cbind(covariates_table.YRI, top.genotypes.YRI)
covariates_table.EUR <- cbind(covariates_table.EUR, top.genotypes.EUR)




    


# RUN DESEQ

# YRI SNVs
Run_DESeq(all_covariates = covariates_table.YRI, 
          SNV_RsID = c('rs2176598'), 
          target_population = c('YRI'), 
          counts_df = counts.YRI, 
          my.output.dir = dir.output.DESeq, 
          my.output.prefix = c('YRI_rs2176598'))  

Run_DESeq(all_covariates = covariates_table.YRI, 
          SNV_RsID = c('rs9271379'), 
          target_population = c('YRI'), 
          counts_df = counts.YRI, 
          my.output.dir = dir.output.DESeq, 
          my.output.prefix = c('YRI_rs9271379'))  

# EUR SNVs
Run_DESeq(all_covariates = covariates_table.EUR, 
          SNV_RsID = c('rs11635336'), 
          target_population = c('EUR'), 
          counts_df = counts.EUR, 
          my.output.dir = dir.output.DESeq, 
          my.output.prefix = c('EUR_rs11635336'))
    
Run_DESeq(all_covariates = covariates_table.EUR, 
          SNV_RsID = c('rs9271894'), 
          target_population = c('EUR'), 
          counts_df = counts.EUR, 
          my.output.dir = dir.output.DESeq, 
          my.output.prefix = c('EUR_rs9271894'))    

Run_DESeq(all_covariates = covariates_table.EUR, 
          SNV_RsID = c('rs1061810'), 
          target_population = c('EUR'), 
          counts_df = counts.EUR, 
          my.output.dir = dir.output.DESeq, 
          my.output.prefix = c('EUR_rs1061810'))    

Run_DESeq(all_covariates = covariates_table.EUR, 
          SNV_RsID = c('rs112581165'), 
          target_population = c('EUR'), 
          counts_df = counts.EUR, 
          my.output.dir = dir.output.DESeq, 
          my.output.prefix = c('EUR_rs112581165'))    

Run_DESeq(all_covariates = covariates_table.EUR, 
          SNV_RsID = c('rs72691418'), 
          target_population = c('EUR'), 
          counts_df = counts.EUR, 
          my.output.dir = dir.output.DESeq, 
          my.output.prefix = c('EUR_rs72691418'))    
    
Run_DESeq(all_covariates = covariates_table.EUR, 
          SNV_RsID = c('rs9270493'), 
          target_population = c('EUR'), 
          counts_df = counts.EUR, 
          my.output.dir = dir.output.DESeq, 
          my.output.prefix = c('EUR_rs9270493'))   





# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/1_DESeq_Genotypes/Session_Info/'
    
    
sink(file = paste(dir.session_info,"Session_Info_GEUVADIS_SNV_DESeq.txt", sep =""))
sessionInfo()
sink()      
    



    
# Clean the environment
rm(list=ls())


