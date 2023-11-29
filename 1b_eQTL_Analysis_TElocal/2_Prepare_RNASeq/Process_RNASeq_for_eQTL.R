# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Process_RNASeq_for_eQTL_functions.R")

# Load libraries
library(edgeR)
library(genefilter)

library(DESeq2) # For VST
library(BiocParallel) # Used by DESeq2 for parallelization
  register(MulticoreParam(6)) # Use six cores
library(limma) # for batch effect removal
library(RNOmni) # for inverse normal transform

  
# Define output directory
dir.output <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/'

    
    
    
# Section 1: FILTER LOW EXPRESSION GENES ------------------------------------------------------------------------------------------------------------------------
 


  
# Load sample meta-data
GEUV_sample_info <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/0_All_Sample_Metadata/E-GEUV-1-unique_and_in_Phase3_snps_TE_SVs.txt", header = T, stringsAsFactors = F, sep = '\t')
rownames(GEUV_sample_info) <- GEUV_sample_info$Comment.ENA_RUN.

# Define directories and lists for RNAseq count tables.
count_Table.dir.YRI <- c("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/1_Read_Counting/TElocal_counts_YRI/")
count_Table.dir.EUR <- c("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/1_Read_Counting/TElocal_counts_EUR/")

cntTable_complete.list.YRI <- list.files(count_Table.dir.YRI, "\\.cntTable$", recursive=FALSE, full.names=TRUE)
cntTable_complete.list.EUR <- list.files(count_Table.dir.EUR, "\\.cntTable$", recursive=FALSE, full.names=TRUE)

# Run function to combine the individual counts tables.
gene_TE_counts_raw_YRI <- aggregate_counts(cntTable_complete.list = cntTable_complete.list.YRI)
gene_TE_counts_raw_EUR <- aggregate_counts(cntTable_complete.list = cntTable_complete.list.EUR)

# Remove genes/TEs with all zeros
gene_TE_counts_raw_YRI <- gene_TE_counts_raw_YRI[which(rowSums(gene_TE_counts_raw_YRI[, -c(1)]) > 0), ]
gene_TE_counts_raw_EUR <- gene_TE_counts_raw_EUR[which(rowSums(gene_TE_counts_raw_EUR[, -c(1)]) > 0), ]

# Shorten sample column names 
colnames(gene_TE_counts_raw_YRI) <- sub("Aligned.sortedByCoord.out.bam*", "", colnames(gene_TE_counts_raw_YRI))
colnames(gene_TE_counts_raw_EUR) <- sub("Aligned.sortedByCoord.out.bam*", "", colnames(gene_TE_counts_raw_EUR))

colnames(gene_TE_counts_raw_YRI) <- sub("*fastp_", "", colnames(gene_TE_counts_raw_YRI))
colnames(gene_TE_counts_raw_EUR) <- sub("*fastp_", "", colnames(gene_TE_counts_raw_EUR))

# change column names to match the names used in the genotyping files.
colnames(gene_TE_counts_raw_YRI) <- GEUV_sample_info[colnames(gene_TE_counts_raw_YRI), 1]
colnames(gene_TE_counts_raw_EUR) <- GEUV_sample_info[colnames(gene_TE_counts_raw_EUR), 1]

# Load repeat subsets
repeats.intergenic.nearby <- read.csv("/Users/juanb/Documents/Bioinformatic_Tools/STAR_genome_indices/TElocal_rmsk_files/GRCh38_GENCODE_rmsk_TE_subsets_vs_GENCODE33/Repeats_intergenic_within_5kb.gtf", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)
repeats.intergenic.distal <- read.csv("/Users/juanb/Documents/Bioinformatic_Tools/STAR_genome_indices/TElocal_rmsk_files/GRCh38_GENCODE_rmsk_TE_subsets_vs_GENCODE33/Repeats_intergenic_farther_than_5kb.gtf", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)
repeats.intronic <- read.csv("/Users/juanb/Documents/Bioinformatic_Tools/STAR_genome_indices/TElocal_rmsk_files/GRCh38_GENCODE_rmsk_TE_subsets_vs_GENCODE33/Repeats_intronic.gtf", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)
repeats.exonic <- read.csv("/Users/juanb/Documents/Bioinformatic_Tools/STAR_genome_indices/TElocal_rmsk_files/GRCh38_GENCODE_rmsk_TE_subsets_vs_GENCODE33/Repeats_exonic.gtf", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)

    # Assign repeat_name to rownames
    rownames(repeats.intergenic.nearby) <- repeats.intergenic.nearby$V9
    rownames(repeats.intergenic.distal) <- repeats.intergenic.distal$V9
    rownames(repeats.intronic) <- repeats.intronic$V9
    rownames(repeats.exonic) <- repeats.exonic$V9
    
# Combine TE loci by subfamily, with no other filters (to compare with TEtranscripts quantification)
gene_TE_counts_raw_YRI.aggregate <- prefilter_TE_aggregation(gene_TE_counts = gene_TE_counts_raw_YRI, 
                                                             aggregate_all = TRUE)

gene_TE_counts_raw_EUR.aggregate <- prefilter_TE_aggregation(gene_TE_counts = gene_TE_counts_raw_EUR, 
                                                             aggregate_all = TRUE)

# Combine TE loci by subfamily, but also splitting depending on location relative to genes (intergenic/exonic/intronic)
gene_TE_counts_raw_YRI.aggregate.split <- prefilter_TE_aggregation(gene_TE_counts = gene_TE_counts_raw_YRI, 
                                                                   aggregate_all = FALSE, 
                                                                   intergenic_nearby_repeats = repeats.intergenic.nearby, 
                                                                   intergenic_distal_repeats = repeats.intergenic.distal, 
                                                                   intronic_repeats = repeats.intronic, 
                                                                   exonic_repeats = repeats.exonic)

gene_TE_counts_raw_EUR.aggregate.split <- prefilter_TE_aggregation(gene_TE_counts = gene_TE_counts_raw_EUR, 
                                                                   aggregate_all = FALSE, 
                                                                   intergenic_nearby_repeats = repeats.intergenic.nearby, 
                                                                   intergenic_distal_repeats = repeats.intergenic.distal, 
                                                                   intronic_repeats = repeats.intronic, 
                                                                   exonic_repeats = repeats.exonic)

# Run function to remove gene version info, combine same genes, and filter lowly expressed genes. *IMPORTANT: PSEUDOAUTOSOMAL GENES (PAR_Y) ARE COMBINED
gene_TE_counts_filtered_YRI.aggregate <- cleanup_and_filter_counts(count_data = gene_TE_counts_raw_YRI.aggregate, filter = TRUE, min_counts_per_sample = 10, fraction_of_samples = 78/86) 
gene_TE_counts_filtered_YRI.aggregate.split <- cleanup_and_filter_counts(count_data = gene_TE_counts_raw_YRI.aggregate.split, filter = TRUE, min_counts_per_sample = 10, fraction_of_samples = 78/86) 

gene_TE_counts_filtered_EUR.aggregate <- cleanup_and_filter_counts(count_data = gene_TE_counts_raw_EUR.aggregate, filter = TRUE, min_counts_per_sample = 10, fraction_of_samples = 323/358) 
gene_TE_counts_filtered_EUR.aggregate.split <- cleanup_and_filter_counts(count_data = gene_TE_counts_raw_EUR.aggregate.split, filter = TRUE, min_counts_per_sample = 10, fraction_of_samples = 323/358) 

    # Save counts files.
    write.table(gene_TE_counts_filtered_YRI.aggregate, file = paste(dir.output, "TElocal_counts_YRI_86_filtered_aggregate", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)
    write.table(gene_TE_counts_filtered_YRI.aggregate.split, file = paste(dir.output, "TElocal_counts_YRI_86_filtered_aggregate_split", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)    
    
    write.table(gene_TE_counts_filtered_EUR.aggregate, file = paste(dir.output, "TElocal_counts_EUR_358_filtered_aggregate", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)
    write.table(gene_TE_counts_filtered_EUR.aggregate.split, file = paste(dir.output, "TElocal_counts_EUR_358_filtered_aggregate_split", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)
    
    
# Remove unneeded files
rm(gene_TE_counts_raw_YRI, gene_TE_counts_raw_EUR, repeats.intergenic.nearby, repeats.intergenic.distal, repeats.intronic, repeats.exonic, gene_TE_counts_raw_YRI.aggregate, gene_TE_counts_raw_EUR.aggregate, gene_TE_counts_raw_YRI.aggregate.split, gene_TE_counts_raw_EUR.aggregate.split)
rm(gene_TE_counts_filtered_YRI.aggregate, gene_TE_counts_filtered_YRI.aggregate.split, gene_TE_counts_filtered_EUR.aggregate, gene_TE_counts_filtered_EUR.aggregate.split)
  

# Section 2: BATCH EFFECT REMOVAL -------------------------------------

      
  
      
# DEFINE KNOWN COVARIATES
      
# Load sample meta data
GEUV_sample_info<- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/0_All_Sample_Metadata/E-GEUV-1-unique_and_in_Phase3_snps_TE_SVs.txt", header = T, stringsAsFactors = F, sep = '\t')
rownames(GEUV_sample_info) <- GEUV_sample_info$Source.Name

# Load genotype PCA covariates
Principal_components <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/Population_Structure_Analysis/COVARIATE_SNV_Population_Structure_PCA.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Load the L1/Alu insertion/deletion covariate table
L1_Alu_insertions_deletions <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/Structural_Variant_Counts/L1_Alu_Insertions_Deletions_Per_Sample.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')
    
# Extract known covariates from the meta data
covariates_table <- GEUV_sample_info[, c('Source.Name','Characteristics.sex.', 'Characteristics.ancestry.category.', 'Performer')]

# Update covariate names
colnames(covariates_table) <- c('id', 'sex', 'ancestry', 'lab') 

# Fill in genotype PCs 1-2 (note: these were calculated separately in each superpopulation)
covariates_table$PC1 <- Principal_components[rownames(covariates_table), 'PC1']
covariates_table$PC2 <- Principal_components[rownames(covariates_table), 'PC2']

# Add the net L1/Alu copy number (ie Alu/L1 insertions - Alu/L1 deletions)
covariates_table$Net_L1_Alu_Copies <- L1_Alu_insertions_deletions[rownames(covariates_table), 'Net_Both']

# Split covariates into YRI and EUR
covariates_table.YRI <- covariates_table[which(covariates_table$ancestry == 'Yoruba'), ]
covariates_table.EUR <- covariates_table[which(covariates_table$ancestry != 'Yoruba'), ]

# Scale and center the net copy number in each population
covariates_table.YRI$Net_L1_Alu_Copies <- scale(covariates_table.YRI$Net_L1_Alu_Copies, center = TRUE, scale = TRUE)
covariates_table.EUR$Net_L1_Alu_Copies <- scale(covariates_table.EUR$Net_L1_Alu_Copies, center = TRUE, scale = TRUE)

    

    
    
    
# VST TRANSFORMATION OF THE COUNT DATA

# Load expression data
gene_TE_counts.YRI <- as.matrix(read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/TElocal_counts_YRI_86_filtered_aggregate_split.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t'))
gene_TE_counts.EUR <- as.matrix(read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/TElocal_counts_EUR_358_filtered_aggregate_split.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t'))
 
# Collect covariates in a dataframe
SampleInfo.YRI <- data.frame(row.names = colnames(gene_TE_counts.YRI),
                             sex = covariates_table.YRI$sex,
                             lab = covariates_table.YRI$lab,
                             PC1 = covariates_table.YRI$PC1,
                             PC2 = covariates_table.YRI$PC2,
                             Net_TE_copies = covariates_table.YRI$Net_L1_Alu_Copies)
                             
SampleInfo.EUR <- data.frame(row.names = colnames(gene_TE_counts.EUR),
                             sex = covariates_table.EUR$sex,
                             ancestry = covariates_table.EUR$ancestry,
                             lab = covariates_table.EUR$lab,
                             PC1 = covariates_table.EUR$PC1,
                             PC2 = covariates_table.EUR$PC2,
                             Net_TE_copies = covariates_table.EUR$Net_L1_Alu_Copies) 
  
# Create DESeq2 object 
dds.YRI <- DESeqDataSetFromMatrix(countData = gene_TE_counts.YRI,
                                  colData = SampleInfo.YRI,
                                  design = ~ sex + lab + PC1 + PC2 + Net_TE_copies)

dds.EUR <- DESeqDataSetFromMatrix(countData = gene_TE_counts.EUR,
                                  colData = SampleInfo.EUR,
                                  design = ~ sex + ancestry + lab + PC1 + PC2 + Net_TE_copies)
  
# Run DESeq2
dds.YRI <- DESeq(dds.YRI, parallel = TRUE)
dds.EUR <- DESeq(dds.EUR, parallel = TRUE)

# Get VST data (includes size factor normalization)
VST.YRI <- assay(varianceStabilizingTransformation(dds.YRI, blind = FALSE))
VST.EUR <- assay(varianceStabilizingTransformation(dds.EUR, blind = FALSE))

    # Save VST data
    write.table(VST.YRI, file = paste(dir.output, 'TElocal_counts_YRI_86_filtered_aggregate_split_VST', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
    write.table(VST.EUR, file = paste(dir.output, 'TElocal_counts_EUR_358_filtered_aggregate_split_VST', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
        
    
    
    
    
    
    
# OBTAIN EBV EXPRESSION AS AN ADDITIONAL COVARIATE
    
# Load VST data
VST.YRI <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/TElocal_counts_YRI_86_filtered_aggregate_split_VST.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')
VST.EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/TElocal_counts_EUR_358_filtered_aggregate_split_VST.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')
    
# Add EBV VST expression to the covariates table    
covariates_table.YRI$EBV_expr_VST <- t(VST.YRI["EBV_gene", rownames(covariates_table.YRI)])
covariates_table.EUR$EBV_expr_VST <- t(VST.EUR["EBV_gene", rownames(covariates_table.EUR)])
    
# Scale and center EBV expression
covariates_table.YRI$EBV_expr_VST <- scale(covariates_table.YRI$EBV_expr_VST, center = TRUE, scale = TRUE)
covariates_table.EUR$EBV_expr_VST <- scale(covariates_table.EUR$EBV_expr_VST, center = TRUE, scale = TRUE)
    
# Remove redundant variables
covariates_table.YRI <- covariates_table.YRI[ , -c(3)] # NOTE: DOUBLE CHECK THE COLUMNS REMOVED. REMOVE ANCESTRY.
covariates_table.EUR <- covariates_table.EUR[, ] # NOTE: DOUBLE CHECK COLUMNS REMOVED.

# Save covariates tables
write.table(covariates_table.YRI, file = paste(dir.output, "COVARIATES_YRI", ".txt", sep =""), sep = "\t" , row.names = T, col.names = NA, quote=F)
write.table(covariates_table.EUR, file = paste(dir.output, "COVARIATES_EUR", ".txt", sep =""), sep = "\t" , row.names = T, col.names = NA, quote=F)

    

    
    
    

# BATCH EFFECT REMOVAL WITH LIMMA

# Load VST expression data
VST.YRI <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/TElocal_counts_YRI_86_filtered_aggregate_split_VST.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')
VST.EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/TElocal_counts_EUR_358_filtered_aggregate_split_VST.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Load the covariates data
covariates_table.YRI  <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/COVARIATES_YRI.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
covariates_table.YRI <- covariates_table.YRI[colnames(VST.YRI), ]      
      
covariates_table.EUR  <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/COVARIATES_EUR.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
covariates_table.EUR <- covariates_table.EUR[colnames(VST.EUR), ]   
    
# Make model matrices
full.model.YRI = model.matrix(~ sex + lab + PC1 + PC2 + Net_L1_Alu_Copies + EBV_expr_VST, data = covariates_table.YRI)
full.model.EUR = model.matrix(~ sex + ancestry + lab + PC1 + PC2 + Net_L1_Alu_Copies + EBV_expr_VST, data = covariates_table.EUR)
           
# From the VST data, regress out lab + ancestry + PCs + sex + TE_Insertions + EBV
batch_corr_YRI <- removeBatchEffect(x = VST.YRI,
                                    batch = NULL,
                                    covariates = full.model.YRI[, c(2:12)], # !!!!!!! Change with SVs
                                    design = full.model.YRI[, c(1), drop = FALSE])

batch_corr_EUR <- removeBatchEffect(x = VST.EUR,
                                    batch = NULL,
                                    covariates = full.model.EUR[, c(2:15)], # !!!!!!! Change with SVs
                                    design = full.model.EUR[, c(1), drop = FALSE])

# Save batch corrected data
write.table(batch_corr_YRI, file = paste(dir.output, 'TElocal_counts_YRI_86_filtered_aggregate_split_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
write.table(batch_corr_EUR, file = paste(dir.output, 'TElocal_counts_EUR_358_filtered_aggregate_split_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)


            

            
            
            
# Section 3: INVERSE NORMAL TRANSFORMATIONS (INT) ---------------------------------------------------------------------------------------------------------------



# Inverse Normal Transformation (INT), so data is normally distributed, as needed for linear regression in eqtl

# Load batch-corrected VST expression data
batch_corr_YRI <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/TElocal_counts_YRI_86_filtered_aggregate_split_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')
batch_corr_EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/TElocal_counts_EUR_358_filtered_aggregate_split_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Run function to calculate INT expression for each gene (THESE VALUES WILL SERVE AS THE INPUT FOR MATRIXEQTL)
INT.vst.regressed.YRI <- Inverse_normal_transform(input.df = batch_corr_YRI)
INT.vst.regressed.EUR <- Inverse_normal_transform(input.df = batch_corr_EUR)

# Run function to save INT matrices, and also generate subsets by gene, L1 subfamilies, and Alu subfamilies
Split_and_save_expression_TElocal_quant(expression_df = INT.vst.regressed.YRI, 
                                        output.dir = dir.output, 
                                        output_file_prefix = 'TElocal_counts_YRI_86_filtered_aggregate_split_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV_INT_', 
                                        gene_prefix = 'ENSG', 
                                        TE_groups = 'Subfamily')

Split_and_save_expression_TElocal_quant(expression_df = INT.vst.regressed.EUR, 
                                        output.dir = dir.output, 
                                        output_file_prefix = 'TElocal_counts_EUR_358_filtered_aggregate_split_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV_INT_', 
                                        gene_prefix = 'ENSG', 
                                        TE_groups = 'Subfamily')


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Session_Info/'
    
sink(file = paste(dir.session_info,"Session_Info_Process_RNASeq_for_eQTL.txt", sep =""))
sessionInfo()
sink()    
    


# Clean the environment
rm(list=ls())


