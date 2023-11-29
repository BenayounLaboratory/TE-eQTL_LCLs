# Set strings as factors
options(stringsAsFactors = F)

# Load libraries  
library(BEDMatrix) # To stream genotypes 
library(biomaRt)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/6_Integrate_cis_trans_eQTLs/Integrate_eQTLs_and_Generate_Plots_Functions.R")     



# Annotate TE trans-eQTLs with gene cis-eQTLs -------------------------------

  
 
# Define general parameters

# Define output dir
output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/4_Integrate_cis_trans_eQTLs/Annotated_and_filtered_eQTLs/'

# Load snp info
all_snp_positions.EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/RESOURCE_SNP_info_EUR.txt", header = T, stringsAsFactors = F, sep = '\t')

    # Assign rownames for easier manipulation
    rownames(all_snp_positions.EUR) <- all_snp_positions.EUR$snpid

# Load INT Gene/TE expression
INT.expression.EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/TElocal_counts_EUR_358_filtered_aggregate_split_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV_INT_All_Genes_TEs.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')





# Annotate intronic L1 eQTLs
  
# Load L1 trans-eQTLs passing FDR thresholding
L1.eQTLs <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Significant_SNVs/L1_intronic_trans_eQTLs_SIGNIFICANT_EUR_2023-11-15.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')

# Load L1 trans-eQTL clumping results
L1.clumps <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Significant_SNVs/L1_intronic_trans_eQTLs_SIGNIFICANT_EUR.clumped", header = TRUE, row.names = NULL, stringsAsFactors = F, sep = "")

# Load cis-eQTLs passing FDR thresholding
cis.eQTLs <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Significant_SNVs/Gene_TElocal_cis_eQTLs_SIGNIFICANT_EUR_2023-11-23.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')

    # Only keep cis-eQTL snps that are trans-eQTLs (use the snps column, a custom snp ID that takes into account duplicate RsIDs)
    cis.eQTLs <- cis.eQTLs[which(cis.eQTLs$snps %in% L1.eQTLs$snps), ]
    
# Link TEs and Genes using the shared snps
combined_cis_trans <- link_cis_and_trans_eQTLs(cis.sig = cis.eQTLs, 
                                               trans.sig = L1.eQTLs, 
                                               gene_TE_INT_expression = INT.expression.EUR)

    # Save data
    write.table(combined_cis_trans, file = paste(output.dir, "EUR_L1_intronic_eQTLs_annotated_0_filters", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
    
# Filter on linear regression FDR
combined_cis_trans_1_filter <- combined_cis_trans[which(combined_cis_trans$Gene_TE_FDR < 0.05), ] 

    # Save data
    write.table(combined_cis_trans_1_filter, file = paste(output.dir, "EUR_L1_intronic_eQTLs_annotated_1_filters_geneTEcorrelation", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
    
# Additionally, filter by whether snps are independent (based on Plink R^2)
combined_cis_trans_2_filter <- combined_cis_trans_1_filter[which(combined_cis_trans_1_filter$RsID %in% L1.clumps$SNP), ] 

    # Save data
    write.table(combined_cis_trans_2_filter, file = paste(output.dir, "EUR_L1_intronic_eQTLs_annotated_2_filters_geneTEcorrelation_LD", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)


    
    
    
# Annotate distal intergenic L1 eQTLs
  
# Load L1 trans-eQTLs passing FDR thresholding
L1.eQTLs <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Significant_SNVs/L1_intergenic_distal_trans_eQTLs_SIGNIFICANT_EUR_2023-11-15.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')

# Load L1 trans-eQTL clumping results
L1.clumps <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Significant_SNVs/L1_intergenic_distal_trans_eQTLs_SIGNIFICANT_EUR.clumped", header = TRUE, row.names = NULL, stringsAsFactors = F, sep = "")

# Load cis-eQTLs passing FDR thresholding
cis.eQTLs <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Significant_SNVs/Gene_TElocal_cis_eQTLs_SIGNIFICANT_EUR_2023-11-23.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')

    # Only keep cis-eQTL snps that are trans-eQTLs (use the snps column, a custom snp ID that takes into account duplicate RsIDs)
    cis.eQTLs <- cis.eQTLs[which(cis.eQTLs$snps %in% L1.eQTLs$snps), ]
    
# Link TEs and Genes using the shared snps
combined_cis_trans <- link_cis_and_trans_eQTLs(cis.sig = cis.eQTLs, 
                                               trans.sig = L1.eQTLs, 
                                               gene_TE_INT_expression = INT.expression.EUR)

    # Save data
    write.table(combined_cis_trans, file = paste(output.dir, "EUR_L1_intergenic_distal_eQTLs_annotated_0_filters", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
    
# Filter on linear regression FDR
combined_cis_trans_1_filter <- combined_cis_trans[which(combined_cis_trans$Gene_TE_FDR < 0.05), ] 

    # Save data
    write.table(combined_cis_trans_1_filter, file = paste(output.dir, "EUR_L1_intergenic_distal_eQTLs_annotated_1_filters_geneTEcorrelation", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
    
# Additionally, filter by whether snps are independent (based on Plink R^2)
combined_cis_trans_2_filter <- combined_cis_trans_1_filter[which(combined_cis_trans_1_filter$RsID %in% L1.clumps$SNP), ] 

    # Save data
    write.table(combined_cis_trans_2_filter, file = paste(output.dir, "EUR_L1_intergenic_distal_eQTLs_annotated_2_filters_geneTEcorrelation_LD", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)

    
    
    
    
# Annotate nearby intergenic L1 eQTLs
  
# Load L1 trans-eQTLs passing FDR thresholding
L1.eQTLs <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Significant_SNVs/L1_intergenic_nearby_trans_eQTLs_SIGNIFICANT_EUR_2023-11-15.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')

# Load L1 trans-eQTL clumping results
L1.clumps <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Significant_SNVs/L1_intergenic_nearby_trans_eQTLs_SIGNIFICANT_EUR.clumped", header = TRUE, row.names = NULL, stringsAsFactors = F, sep = "")

# Load cis-eQTLs passing FDR thresholding
cis.eQTLs <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Significant_SNVs/Gene_TElocal_cis_eQTLs_SIGNIFICANT_EUR_2023-11-23.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')

    # Only keep cis-eQTL snps that are trans-eQTLs (use the snps column, a custom snp ID that takes into account duplicate RsIDs)
    cis.eQTLs <- cis.eQTLs[which(cis.eQTLs$snps %in% L1.eQTLs$snps), ]
    
# Link TEs and Genes using the shared snps
combined_cis_trans <- link_cis_and_trans_eQTLs(cis.sig = cis.eQTLs, 
                                               trans.sig = L1.eQTLs, 
                                               gene_TE_INT_expression = INT.expression.EUR)

    # Save data
    write.table(combined_cis_trans, file = paste(output.dir, "EUR_L1_intergenic_nearby_eQTLs_annotated_0_filters", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
    
# Filter on linear regression FDR
combined_cis_trans_1_filter <- combined_cis_trans[which(combined_cis_trans$Gene_TE_FDR < 0.05), ] 

    # Save data
    write.table(combined_cis_trans_1_filter, file = paste(output.dir, "EUR_L1_intergenic_nearby_eQTLs_annotated_1_filters_geneTEcorrelation", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
    
# Additionally, filter by whether snps are independent (based on Plink R^2)
combined_cis_trans_2_filter <- combined_cis_trans_1_filter[which(combined_cis_trans_1_filter$RsID %in% L1.clumps$SNP), ] 

    # Save data
    write.table(combined_cis_trans_2_filter, file = paste(output.dir, "EUR_L1_intergenic_nearby_eQTLs_annotated_2_filters_geneTEcorrelation_LD", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)

    
    
    
    
# Annotate exonic L1 eQTLs
  
# Load L1 trans-eQTLs passing FDR thresholding
L1.eQTLs <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Significant_SNVs/L1_exonic_trans_eQTLs_SIGNIFICANT_EUR_2023-11-15.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')

# Load L1 trans-eQTL clumping results
L1.clumps <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Significant_SNVs/L1_exonic_trans_eQTLs_SIGNIFICANT_EUR.clumped", header = TRUE, row.names = NULL, stringsAsFactors = F, sep = "")

# Load cis-eQTLs passing FDR thresholding
cis.eQTLs <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/3_Scan_for_eQTLs/Significant_SNVs/Gene_TElocal_cis_eQTLs_SIGNIFICANT_EUR_2023-11-23.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')

    # Only keep cis-eQTL snps that are trans-eQTLs (use the snps column, a custom snp ID that takes into account duplicate RsIDs)
    cis.eQTLs <- cis.eQTLs[which(cis.eQTLs$snps %in% L1.eQTLs$snps), ]
    
# Link TEs and Genes using the shared snps
combined_cis_trans <- link_cis_and_trans_eQTLs(cis.sig = cis.eQTLs, 
                                               trans.sig = L1.eQTLs, 
                                               gene_TE_INT_expression = INT.expression.EUR)

    # Save data
    write.table(combined_cis_trans, file = paste(output.dir, "EUR_L1_exonic_eQTLs_annotated_0_filters", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
    
# Filter on linear regression FDR
combined_cis_trans_1_filter <- combined_cis_trans[which(combined_cis_trans$Gene_TE_FDR < 0.05), ] 

    # Save data
    write.table(combined_cis_trans_1_filter, file = paste(output.dir, "EUR_L1_exonic_eQTLs_annotated_1_filters_geneTEcorrelation", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
    
# Additionally, filter by whether snps are independent (based on Plink R^2)
combined_cis_trans_2_filter <- combined_cis_trans_1_filter[which(combined_cis_trans_1_filter$RsID %in% L1.clumps$SNP), ] 

    # Save data
    write.table(combined_cis_trans_2_filter, file = paste(output.dir, "EUR_L1_exonic_eQTLs_annotated_2_filters_geneTEcorrelation_LD", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)

    
          


      
  

      
      
      
        
# Plot snp-gene-TE trio relationships (INTERGENIC DISTAL) -------------------------------        
        
        
      
    
    
# Define output dir
output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/4_Integrate_cis_trans_eQTLs/Genotype_TE_Gene_Plots/'

# Load file with snp info
all_snp_positions <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/RESOURCE_SNP_info_EUR.txt", header = T, stringsAsFactors = F, sep = '\t')

    # Update rownames with snpid for easier access
    rownames(all_snp_positions) <- all_snp_positions$snpid

# Specify the genotype BED file
plink_file_path <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink.bed'

    # Make R object to stream snp genotypes into memory (rows are samples and columns are snps). # NOTE: If used, MT snps will need to be assigned a snpid, since they dont have rsid.
    binary_genotypes <- BEDMatrix(path = plink_file_path, simple_names = T)

# Load Gene/TE expression
Gene_expression <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/TElocal_counts_EUR_358_filtered_aggregate_split_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV_INT_All_Genes_TEs.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Load combined, cis-trans QTL table
combined_cis_trans <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/4_Integrate_cis_trans_eQTLs/Annotated_and_filtered_eQTLs/EUR_L1_intergenic_distal_eQTLs_annotated_2_filters_geneTEcorrelation_LD.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')

    # For entries with no symbol, use the EnsemblID
    combined_cis_trans[which(combined_cis_trans$symbol == ""), 'symbol'] <- combined_cis_trans[which(combined_cis_trans$symbol == ""), 'gene']

# Generate plots as above, but ONLY for associations with gene symbols
eQTL_to_barplots_scatterplot(output.dir = output.dir,
                             output.filename = 'EUR_L1_Intergenic_Distal_eQTLs_Primary_Candidates_SNP1_Symbols_only_',
                             plot_rows = 3,
                             plot_columns = 3,
                             indices_to_plot = c(1:1),
                             combined_cis_trans_QTL_df = combined_cis_trans,
                             all_snp_positions = all_snp_positions,
                             gene_expression_df = Gene_expression,
                             streamed_genotypes_matrix = binary_genotypes)
    
        

# Clean the environment
rm(list=ls())       


