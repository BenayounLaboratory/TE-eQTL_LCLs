# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/Run_MatrixeQTL_Functions.R") 

# Load libraries
library(MatrixEQTL) # package that will be carrying out the eQTL analysis


  
#### EQTLS WITH EUR (BOTH SEXES) POPULATIONS

# Define genotype file (complete, not pruned) and load it (to avoid reloading during each eQTL analysis)
SNP_file_name = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink.012"
snps <- load_eQTL_genotypes(expression_file = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/TElocal_counts_EUR_358_filtered_aggregate_split_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV_INT_All_Genes_TEs.txt",
                            SNP_file = SNP_file_name)  

# intergenic (near genes) L1 trans-eQTLs
Single_Matrix_eQTL(is.perm = FALSE,
                   output_file_cis = NULL,
                   snps_pos = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/RESOURCE_SNP_info_EUR.txt",
                   gene_pos = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/MatrixeQTL_Additional_Input_Files/ENSEMBL_GENE_POSITIONS_EUR.txt",
                   pvOutputThreshold_cis = 0,
                   cis_Dist = 1e6,
                   output_file_tra = NULL,
                   pvOutputThreshold_tra = 1e-2,
                   expression_file = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/TElocal_counts_EUR_358_filtered_aggregate_split_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV_INT_Only_L1_intergenic_near.txt",
                   gene = NULL,
                   SNP_file = NULL,
                   snps = snps,
                   covariates_file = NULL,
                   R_obj = F, 
                   output.dir = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Raw_eQTL_Results/',
                   R_object_filename = 'L1_intergenic_near_genes_subfamily_trans_EUR')

# intergenic (distal to genes) L1 trans-eQTLs
Single_Matrix_eQTL(is.perm = FALSE,
                   output_file_cis = NULL,
                   snps_pos = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/RESOURCE_SNP_info_EUR.txt",
                   gene_pos = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/MatrixeQTL_Additional_Input_Files/ENSEMBL_GENE_POSITIONS_EUR.txt",
                   pvOutputThreshold_cis = 0,
                   cis_Dist = 1e6,
                   output_file_tra = NULL,
                   pvOutputThreshold_tra = 1e-2,
                   expression_file = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/TElocal_counts_EUR_358_filtered_aggregate_split_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV_INT_Only_L1_intergenic_distal.txt",
                   gene = NULL,
                   SNP_file = NULL,
                   snps = snps,
                   covariates_file = NULL,
                   R_obj = F, 
                   output.dir = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Raw_eQTL_Results/',
                   R_object_filename = 'L1_intergenic_distal_genes_subfamily_trans_EUR')

# intronic L1 trans-eQTLs
Single_Matrix_eQTL(is.perm = FALSE,
                   output_file_cis = NULL,
                   snps_pos = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/RESOURCE_SNP_info_EUR.txt",
                   gene_pos = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/MatrixeQTL_Additional_Input_Files/ENSEMBL_GENE_POSITIONS_EUR.txt",
                   pvOutputThreshold_cis = 0,
                   cis_Dist = 1e6,
                   output_file_tra = NULL,
                   pvOutputThreshold_tra = 1e-2,
                   expression_file = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/TElocal_counts_EUR_358_filtered_aggregate_split_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV_INT_Only_L1_intronic.txt",
                   gene = NULL,
                   SNP_file = NULL,
                   snps = snps,
                   covariates_file = NULL,
                   R_obj = F, 
                   output.dir = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Raw_eQTL_Results/',
                   R_object_filename = 'L1_intronic_subfamily_trans_EUR')

# exonic L1 trans-eQTLs
Single_Matrix_eQTL(is.perm = FALSE,
                   output_file_cis = NULL,
                   snps_pos = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/RESOURCE_SNP_info_EUR.txt",
                   gene_pos = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/MatrixeQTL_Additional_Input_Files/ENSEMBL_GENE_POSITIONS_EUR.txt",
                   pvOutputThreshold_cis = 0,
                   cis_Dist = 1e6,
                   output_file_tra = NULL,
                   pvOutputThreshold_tra = 1e-2,
                   expression_file = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/TElocal_counts_EUR_358_filtered_aggregate_split_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV_INT_Only_L1_exonic.txt",
                   gene = NULL,
                   SNP_file = NULL,
                   snps = snps,
                   covariates_file = NULL,
                   R_obj = F, 
                   output.dir = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Raw_eQTL_Results/',
                   R_object_filename = 'L1_exonic_subfamily_trans_EUR')

# Gene cis-eQTLs
Single_Matrix_eQTL(is.perm = FALSE,
                     output_file_cis = NULL,
                     snps_pos = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/RESOURCE_SNP_info_EUR.txt",
                     gene_pos = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/MatrixeQTL_Additional_Input_Files/ENSEMBL_GENE_POSITIONS_EUR.txt",
                     pvOutputThreshold_cis = 1e-1,
                     cis_Dist = 1e6,
                     output_file_tra = NULL,
                     pvOutputThreshold_tra = 0,
                     expression_file = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/TElocal_counts_EUR_358_filtered_aggregate_split_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV_INT_Only_Genes.txt",
                     gene = NULL,
                     SNP_file = NULL,
                     snps = snps,
                     covariates_file = NULL,
                     R_obj = F, 
                     output.dir = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Raw_eQTL_Results/',
                     R_object_filename = 'Gene_cis_EUR_TElocal')




# # PERMUTATIONS: L1 trans-eQTLs
# Permuted_eQTLs(permutations_file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/MatrixeQTL_Additional_Input_Files/Input_Permutations_Both_Sexes_EUR.txt',
#                specific_permutations = 1:20,
#                output_file_cis = NULL,
#                snps_pos = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/RESOURCE_SNP_info_EUR.txt",
#                gene_pos = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/MatrixeQTL_Additional_Input_Files/ENSEMBL_GENE_POSITIONS_EUR.txt",
#                pvOutputThreshold_cis = 0,
#                cis_Dist = 1e6,
#                output_file_tra = NULL,
#                pvOutputThreshold_tra = 1e-2,
#                expression_file = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_eQTL/All_counts_EUR_358_filtered_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV_INT_Only_L1.txt",
#                snps = snps,
#                covariates_file = NULL,
#                output.dir = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/Raw_eQTL_Results/L1_subfamily_trans_EUR_Permutations/',
#                R_object_filename = 'L1_subfamily_trans')



  



    
    
## Clean the environment
rm(list=ls())    

