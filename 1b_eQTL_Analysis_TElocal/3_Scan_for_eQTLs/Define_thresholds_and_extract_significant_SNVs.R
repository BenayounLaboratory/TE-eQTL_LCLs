# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/Define_thresholds_and_extract_significant_SNVs_functions.R") 

# Define output directory for significant snp files
output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Significant_SNVs/'



# Determine Empirical FDR/pvalue thresholds --------------------------------------------



# # Define input parameters
# 
# # Define the target empirical FDR threshold
# empirical.FDR.threshold <- 0.05
# 
# # Make a dataframe to hold the pvalue thresholds (and other stats) for the desired empirical FDR threshold
# all.empirical.data <- as.data.frame(matrix(nrow = 2, ncol = 4, NA))
# 
# # Assign colnames
# colnames(all.empirical.data) <- c('P_threshold', 'Empirical_FDR_at_P', 'Avg_Num_Permutation_Points_Below_P', 'Avg_Num_Real_Points_Below_P')
# 
# # Assign rownames
# rownames(all.empirical.data) <- c('EUR_L1_trans', 'EUR_Gene_cis')
# 
# 
# 
# 
# 
# # Calculate empirical FDR thresholds 
#   
# # EUR L1 trans-eQTL
# empirical_summary_stats <- empiricalFDR_to_pval(data.file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/Raw_eQTL_Results/L1_subfamily_trans_EUR.rds', 
#                                                 null_data.dir = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/Raw_eQTL_Results/L1_subfamily_trans_EUR_Permutations/', 
#                                                 eQTL_type = 'trans', 
#                                                 arbitrary_pval = 1e-7, 
#                                                 empirical_FDR_target = empirical.FDR.threshold)
# 
#     # Print out the results
#     empirical_summary_stats # 3.442833e-08 4.878731e-02 2.615000e+01 5.360000e+02
#     
#     # Add results to the df
#     if (!is.null(empirical_summary_stats)) {
#         all.empirical.data['EUR_L1_trans', ] <- empirical_summary_stats
#     }
#       
# 
#     
#     
# 
# 
# 
# 
# # Save the results
# write.table(all.empirical.data, file = paste(output.dir, "EMPIRICAL_PVALUE_AND_FDR_THRESHOLDS", '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)

    
    
        


# Extract Significant SNPs --------------------------------------------



          
          
# Define various input parameters

# Load table with EUR snp details
all_snp_positions.EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/RESOURCE_SNP_info_EUR.txt", header = T, stringsAsFactors = F, sep = '\t')

# Load table with empirical pval/FDR thresholds
#empirical.thresholds <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/Significant_SNVs/EMPIRICAL_PVALUE_AND_FDR_THRESHOLDS.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Define the desired BH FDR threshold 
BH.FDR.threshold <- 0.05

  


  
# EXTRACT SIGNIFICANT SNPS : EUR intergenic distal L1s

# Read eQTL results. Only keep stats.
eQTL.stats <- readRDS('/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Raw_eQTL_Results/L1_intergenic_distal_genes_subfamily_trans_EUR.rds')
eQTL.stats <- eQTL.stats$all$eqtls

# Extract pvalue corresponding to the desired BH FDR
BH.pval_threshold <- max(eQTL.stats[eQTL.stats$FDR < BH.FDR.threshold, 'pvalue']) 
BH.pval_threshold # 2.724676e-11

# Define empirical FDR pvalue
#empirical.pval_threshold <- empirical.thresholds[c('EUR_L1_trans'), c('P_threshold')]
  
# Run function to extract significant snps
Extract_sig_snps_from_matrixEQTL(output_location = output.dir, 
                                 eQTL_stats_df = eQTL.stats, 
                                 vector_of_significance_thresholds = c(BH.pval_threshold), 
                                 snp_info_df = all_snp_positions.EUR,
                                 output.file.label = c('L1_intergenic_distal_trans_eQTLs_SIGNIFICANT_EUR_'),
                                 save_rsid_list = c('no'))





# EXTRACT SIGNIFICANT SNPS : EUR intergenic nearby L1s

# Read eQTL results. Only keep stats.
eQTL.stats <- readRDS('/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Raw_eQTL_Results/L1_intergenic_near_genes_subfamily_trans_EUR.rds')
eQTL.stats <- eQTL.stats$all$eqtls

# Extract pvalue corresponding to the desired BH FDR
BH.pval_threshold <- max(eQTL.stats[eQTL.stats$FDR < BH.FDR.threshold, 'pvalue']) 
BH.pval_threshold # 4.704457e-07

# Define empirical FDR pvalue
#empirical.pval_threshold <- empirical.thresholds[c('EUR_L1_trans'), c('P_threshold')]
  
# Run function to extract significant snps
Extract_sig_snps_from_matrixEQTL(output_location = output.dir, 
                                 eQTL_stats_df = eQTL.stats, 
                                 vector_of_significance_thresholds = c(BH.pval_threshold), 
                                 snp_info_df = all_snp_positions.EUR,
                                 output.file.label = c('L1_intergenic_nearby_trans_eQTLs_SIGNIFICANT_EUR_'),
                                 save_rsid_list = c('no'))





# EXTRACT SIGNIFICANT SNPS : EUR intronic L1s

# Read eQTL results. Only keep stats.
eQTL.stats <- readRDS('/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Raw_eQTL_Results/L1_intronic_subfamily_trans_EUR.rds')
eQTL.stats <- eQTL.stats$all$eqtls

# Extract pvalue corresponding to the desired BH FDR
BH.pval_threshold <- max(eQTL.stats[eQTL.stats$FDR < BH.FDR.threshold, 'pvalue']) 
BH.pval_threshold # 2.686097e-08

# Define empirical FDR pvalue
#empirical.pval_threshold <- empirical.thresholds[c('EUR_L1_trans'), c('P_threshold')]
  
# Run function to extract significant snps
Extract_sig_snps_from_matrixEQTL(output_location = output.dir, 
                                 eQTL_stats_df = eQTL.stats, 
                                 vector_of_significance_thresholds = c(BH.pval_threshold), 
                                 snp_info_df = all_snp_positions.EUR,
                                 output.file.label = c('L1_intronic_trans_eQTLs_SIGNIFICANT_EUR_'),
                                 save_rsid_list = c('no'))





# EXTRACT SIGNIFICANT SNPS : EUR exonic L1s

# Read eQTL results. Only keep stats.
eQTL.stats <- readRDS('/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Raw_eQTL_Results/L1_exonic_subfamily_trans_EUR.rds')
eQTL.stats <- eQTL.stats$all$eqtls

# Extract pvalue corresponding to the desired BH FDR
BH.pval_threshold <- max(eQTL.stats[eQTL.stats$FDR < BH.FDR.threshold, 'pvalue']) 
BH.pval_threshold # 1.382019e-08

# Define empirical FDR pvalue
#empirical.pval_threshold <- empirical.thresholds[c('EUR_L1_trans'), c('P_threshold')]
  
# Run function to extract significant snps
Extract_sig_snps_from_matrixEQTL(output_location = output.dir, 
                                 eQTL_stats_df = eQTL.stats, 
                                 vector_of_significance_thresholds = c(BH.pval_threshold), 
                                 snp_info_df = all_snp_positions.EUR,
                                 output.file.label = c('L1_exonic_trans_eQTLs_SIGNIFICANT_EUR_'),
                                 save_rsid_list = c('no'))
            




# EXTRACT SIGNIFICANT SNPS : EUR Gene cis eQTL

# Read external files. Only keep stats.
eQTL.stats <- readRDS('/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Raw_eQTL_Results/Gene_cis_EUR_TElocal.rds')
eQTL.stats <- eQTL.stats$cis$eqtls
  
# Extract pvalue corresponding to user-defined FDR
BH.pval_threshold <- max(eQTL.stats[eQTL.stats$FDR < BH.FDR.threshold, 'pvalue']) 
BH.pval_threshold # 0.000474156

# Define empirical FDR pvalue
# empirical.pval_threshold <- empirical.thresholds[c('EUR_Gene_cis'), c('P_threshold')]

# Run function to extract significant snps
Extract_sig_snps_from_matrixEQTL(output_location = output.dir, 
                                 eQTL_stats_df = eQTL.stats, 
                                 vector_of_significance_thresholds = c(BH.pval_threshold), 
                                 snp_info_df = all_snp_positions.EUR,
                                 output.file.label = c('Gene_TElocal_cis_eQTLs_SIGNIFICANT_EUR_'),
                                 save_rsid_list = c('no'))


  


# Clean the environment
rm(list=ls())       



