# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/Prepare_Remaining_Input_Files_Functions.R") 

# Load libraries
library(arrangements) # to generate permutations

# Specify output directory
dir.output <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/MatrixeQTL_Additional_Input_Files/'

  




## Define ENSEMBL gene positions using GENCODE GTF. 

# Load expressed genes
expression.YRI <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_eQTL/All_counts_YRI_86_filtered_VST.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
expression.EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_eQTL/All_counts_EUR_358_filtered_VST.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Make a copy of the Gencode v30 GTF and delete the header

# In Terminal, make a bed file containing only "gene" entries
# awk '$3 == "gene"' gencode_copy.gtf > gencode_33_genes.bed 

# Step 3: Load the gene bed file
gene_info <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/MatrixeQTL_Additional_Input_Files/gencode_33_genes.bed", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)

# Remove extra info attached to the ENSEMBLE ID
gene_info$V9 <- sub("gene_id *", "", gene_info$V9) # Remove the label before the id
gene_info$V9 <- sub("\\..*", "", gene_info$V9, fixed=FALSE) # Remove everything after the id. NOTE: PAR GENE LABEL REMOVED

# Change chromosome labels (remove 'chr')
gene_info$V1 <- gsub('chr', '', gene_info$V1)

# Keep genes on autosomes
gene_info <- gene_info[which(gene_info$V1 %in% c(1:22)), ]

# Check for gene duplicates
sum(duplicated(gene_info$V9)) == 0

# Keep useful columns
gene_info <- gene_info[, c('V9', 'V1', 'V4', 'V5')]
colnames(gene_info) <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position")
    
# Subset info for expressed genes
    
# Make df for each population
gene_info.YRI <- gene_info
gene_info.EUR <- gene_info

# Check which genes are expressed ** AND will be analyzed (the analysis will not include genes on non-autosomal chromosomes)
to_keep.YRI <- which(gene_info.YRI$ensembl_gene_id %in% rownames(expression.YRI))
to_keep.EUR <- which(gene_info.EUR$ensembl_gene_id %in% rownames(expression.EUR))

# Subset the info for the final genes of interest
gene_info.YRI <- gene_info.YRI[to_keep.YRI, ]
gene_info.EUR <- gene_info.EUR[to_keep.EUR, ]

# Save gene position file
write.table(gene_info.YRI, file = paste(dir.output, "ENSEMBL_GENE_POSITIONS_YRI", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
write.table(gene_info.EUR, file = paste(dir.output, "ENSEMBL_GENE_POSITIONS_EUR", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)

  





# GENERATE SAMPLE PERMUTATION FILES (for empirical pvalue thresholding in eQTL analysis)
    
# Permute samples names (for eQTL null distribution)
Both_sexes_YRI_sample_permutations <- generate_permutations(vector_to_permute = colnames(expression.YRI),
                                                            unique_only = T,
                                                            perm_length = ncol(expression.YRI),
                                                            max_permutations = 1e5,
                                                            perm_layout = 'column',
                                                            perms_to_keep = 1000,
                                                            output_file = paste(dir.output, "Input_Permutations_Both_Sexes_YRI", '.txt', sep =""))

Both_sexes_EUR_sample_permutations <- generate_permutations(vector_to_permute = colnames(expression.EUR),
                                                            unique_only = T,
                                                            perm_length = ncol(expression.EUR),
                                                            max_permutations = 1e5,
                                                            perm_layout = 'column',
                                                            perms_to_keep = 1000,
                                                            output_file = paste(dir.output, "Input_Permutations_Both_Sexes_EUR", '.txt', sep =""))



    
## Clean the environment
rm(list=ls())
      






