# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with or needed by this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/7_Replication_of_eQTLs_in_YRI/Check_Replication_of_eQTLs_Functions.R") 
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/6_Integrate_cis_trans_eQTLs/Integrate_eQTLs_and_Generate_Plots_Functions.R") 
      
# load libraries
library(BEDMatrix) # needed to stream genotypes
      
      
    
  
# Define universal variables

# Define the YRI gene/TE INT expression file path
Gene_expression.YRI.path <- "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_eQTL/All_counts_YRI_86_filtered_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV_INT_All_Genes_TEs.txt"

# Load Gene/TE INT expression
Gene_expression.YRI <- read.csv(Gene_expression.YRI.path, header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Define the YRI genotype file path
YRI.genotype.bed.path <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_YRI_Genotypes/Filtered_YRI_BOTH_SEXES_chr1-22.recode.plink.bed'

    # Make R object to stream snp genotypes into memory (rows are samples and columns are snps). # NOTE: If used, MT snps will need to be assigned a snpid, since they dont have rsid.
    binary_genotypes <- BEDMatrix(path = YRI.genotype.bed.path, simple_names = T)
    
# Load snp details file
all_snp_positions.YRI <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_YRI_Genotypes/RESOURCE_SNP_info_YRI.txt", header = T, stringsAsFactors = F, sep = '\t')

    # Remove duplicate RsIDs.There should only be one. Keep only the first, and use the snpid as a placeholder for the below analyses. 
    all_snp_positions.YRI <- all_snp_positions.YRI[!duplicated(all_snp_positions.YRI$RsID), ]

    # Assign RsID to rownames
    rownames(all_snp_positions.YRI) <- all_snp_positions.YRI$RsID
    
    
    


# SCAN FOR EUR EQTL REPLICATION IN YRI ----------------------------



# Define output directory 
output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/7_Replication_of_eQTLs_in_YRI/Scan_for_eQTLs/'
  
  
    
        
# Run EUR L1 trans-eQTLs
  
# Load EUR significant snps
sig.eQTLs <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/Significant_SNVs/L1_trans_eQTLs_SIGNIFICANT_EUR_2023-05-22.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')

# Define indices of sig. EUR snps that are present in the YRI population
EUR.shared_sig_snp_indices <- which(sig.eQTLs$RsID %in% all_snp_positions.YRI$RsID)

# Keep EUR stats for YRI-shared snps
sig.eQTLs <- sig.eQTLs[EUR.shared_sig_snp_indices, ]

# Determine sig of EUR snps in YRI population
EUR_snps_in_YRI <- check_if_QTLs_replicate_function(snp_info_df = all_snp_positions.YRI,
                                                    genotype_BED_file_path = YRI.genotype.bed.path,
                                                    Expression_file_path = Gene_expression.YRI.path,
                                                    sig_QTL_df = sig.eQTLs)

# Sort results by pval
EUR_snps_in_YRI <- EUR_snps_in_YRI[order(EUR_snps_in_YRI$pvalue, decreasing = FALSE), ]

    # Save unfiltered results
    write.table(EUR_snps_in_YRI, file = paste(output.dir, "L1_subfamily_trans_eQTLs_YRI_ALL_", Sys.Date(), '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)

# Filter by BH FDR
EUR_snps_in_YRI <- EUR_snps_in_YRI[which(EUR_snps_in_YRI$FDR < 0.05), ]

# Save data
write.table(EUR_snps_in_YRI, file = paste(output.dir, "L1_subfamily_trans_eQTLs_YRI_FDR5_", Sys.Date(), '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)

    

    
# Run EUR gene cis-eQTLs
  
# Load EUR significant cis-eQTLs LINKED TO A TE SNV
sig.eQTLs <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/6_Integrate_cis_trans_eQTLs/Annotated_and_filtered_eQTLs/EUR_L1_eQTLs_annotated_0_filters_2023-05-22.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')

    # Remove rows with NA genes
    sig.eQTLs <- na.omit(sig.eQTLs)
    
    # Only keep QTLs with genes that are expressed in YRI
    sig.eQTLs <- sig.eQTLs[which(sig.eQTLs$gene %in% rownames(Gene_expression.YRI)), ]

# Define indices of sig. EUR snps that are present in the YRI population
EUR.shared_sig_snp_indices <- which(sig.eQTLs$RsID %in% all_snp_positions.YRI$RsID)

# Keep stats for YRI shared snps
sig.eQTLs <- sig.eQTLs[EUR.shared_sig_snp_indices, ]

# Determine sig of EUR snps in YRI
EUR_snps_in_YRI <- check_if_QTLs_replicate_function(snp_info_df = all_snp_positions.YRI,
                                                    genotype_BED_file_path = YRI.genotype.bed.path,
                                                    Expression_file_path = Gene_expression.YRI.path,
                                                    sig_QTL_df = sig.eQTLs)

# Sort by pval
EUR_snps_in_YRI <- EUR_snps_in_YRI[order(EUR_snps_in_YRI$pvalue, decreasing = FALSE), ]

# Save unfiltered data
write.table(EUR_snps_in_YRI, file = paste(output.dir, "Gene_cis_eQTLs_YRI_ALL_", Sys.Date(), '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)

# Filter by BH FDR
EUR_snps_in_YRI <- EUR_snps_in_YRI[which(EUR_snps_in_YRI$FDR < 0.05), ]

# Save data
write.table(EUR_snps_in_YRI, file = paste(output.dir, "Gene_cis_eQTLs_YRI_FDR5_", Sys.Date(), '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
   
  

    
# # Run EUR Alu trans-eQTL (DID NOT MAKE IT TO THE FINAL ANALYSIS)
# 
# # Load EUR significant snps
# sig.eQTLs <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/Significant_SNVs/Alu_trans_eQTLs_SIGNIFICANT_EUR_2022-10-26.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')
# 
# # Define indices of sig. EUR snps that are present in the YRI population
# EUR.shared_sig_snp_indices <- which(sig.eQTLs$RsID %in% all_snp_positions.YRI$RsID)
# 
# # Keep stats for YRI shared snps
# sig.eQTLs <- sig.eQTLs[EUR.shared_sig_snp_indices, ]
# 
# # Determine sig of EUR snps in YRI
# EUR_snps_in_YRI <- check_if_QTLs_replicate_function(snp_info_df = all_snp_positions.YRI,
#                                                     genotype_BED_file_path = YRI.genotype.bed.path,
#                                                     Expression_file_path = Gene_expression.YRI.path,
#                                                     sig_QTL_df = sig.eQTLs)
# 
# # Sort by pval
# EUR_snps_in_YRI <- EUR_snps_in_YRI[order(EUR_snps_in_YRI$pvalue, decreasing = FALSE), ]
# 
# # Filter by BH FDR
# EUR_snps_in_YRI <- EUR_snps_in_YRI[which(EUR_snps_in_YRI$FDR < 0.05), ]
# 
# # Save data
# write.table(EUR_snps_in_YRI, file = paste(output.dir, "Alu_subfamily_trans_eQTLs_YRI_FDR5_", Sys.Date(), '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)


        
     
# REMOVE INDEX SNPS IN LD WITH LEAD SNPS, USING PLINK        
# Run the 'Clump_Significant_SNVs.sh' script
        
        
        
        
        
# INTEGRATE CIS AND TRANS QTLs AND PLOT ----------------------------



        
# INTEGRATE RESULTS

# Define output dir
output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/7_Replication_of_eQTLs_in_YRI/Annotated_and_filtered_eQTLs/'
  
# Load YRI gene cis-eQTLs
cis.eQTLs <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/7_Replication_of_eQTLs_in_YRI/Scan_for_eQTLs/Gene_cis_eQTLs_YRI_FDR5_2023-05-23.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')

# Load YRI L1 trans-eQTLs
L1.eQTLs <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/7_Replication_of_eQTLs_in_YRI/Scan_for_eQTLs/L1_subfamily_trans_eQTLs_YRI_FDR5_2023-05-23.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')
        
# Load clumping results
L1.clumps <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/7_Replication_of_eQTLs_in_YRI/Scan_for_eQTLs/L1_subfamily_trans_eQTLs_YRI_FDR5.clumped", header = TRUE, row.names = NULL, stringsAsFactors = F, sep = "")
       
# Link TEs and Genes using the shared snps
combined_cis_trans <- link_cis_and_trans_eQTLs(cis.sig = cis.eQTLs, trans.sig = L1.eQTLs, gene_TE_INT_expression = Gene_expression.YRI)

    # Save data
    write.table(combined_cis_trans, file = paste(output.dir, "YRI_L1_eQTLs_annotated_0_filters_", Sys.Date(), '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
    
# Filter on linear regression FDR
combined_cis_trans_1_filter <- combined_cis_trans[which(combined_cis_trans$Gene_TE_FDR < 0.05), ] 

    # Save data
    write.table(combined_cis_trans_1_filter, file = paste(output.dir, "YRI_L1_eQTLs_annotated_1_filters_geneTEcorrelation_", Sys.Date(), '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
  
# Subset SNVs that are the most significant in a linkage disequilibrium block
combined_cis_trans_2_filter <- combined_cis_trans_1_filter[which(combined_cis_trans_1_filter$RsID %in% L1.clumps$SNP), ] 

    # Save data
    write.table(combined_cis_trans_2_filter, file = paste(output.dir, "YRI_L1_eQTLs_annotated_2_filters_geneTEcorrelation_LD_", Sys.Date(), '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
  
      

        
        
# PLOT GENOTYPE-GENE-TE RELATIONSHIPS        
  
# Define output dir
output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/7_Replication_of_eQTLs_in_YRI/Genotype_TE_Gene_Plots/'

# Update the snp details file by assigning my custom snpid to the rownames
rownames(all_snp_positions.YRI) <- all_snp_positions.YRI$snpid      
    
# Load the integrated, cis-trans eQTL table
combined_cis_trans <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/7_Replication_of_eQTLs_in_YRI/Annotated_and_filtered_eQTLs/YRI_L1_eQTLs_annotated_2_filters_geneTEcorrelation_LD_2023-05-23.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')

    # For entries with no symbol, use the EnsemblID
    combined_cis_trans[which(combined_cis_trans$symbol == ""), 'symbol'] <- combined_cis_trans[which(combined_cis_trans$symbol == ""), 'gene']
    
# Make plots for all associations that pass linear regression and linkage disequilibrium filtering
eQTL_to_barplots_scatterplot(output.dir = output.dir, 
                             output.filename = 'YRI_L1_eQTLs_Primary_Candidates_SNP1_', 
                             plot_rows = 3, 
                             plot_columns = 3, 
                             indices_to_plot = c(1:3), 
                             combined_cis_trans_QTL_df = combined_cis_trans, 
                             all_snp_positions = all_snp_positions.YRI, 
                             gene_expression_df = Gene_expression.YRI, 
                             streamed_genotypes_matrix = binary_genotypes)

eQTL_to_barplots_scatterplot(output.dir = output.dir, 
                             output.filename = 'YRI_L1_eQTLs_Primary_Candidates_SNP2_', 
                             plot_rows = 3, 
                             plot_columns = 3, 
                             indices_to_plot = c(4:4), 
                             combined_cis_trans_QTL_df = combined_cis_trans, 
                             all_snp_positions = all_snp_positions.YRI, 
                             gene_expression_df = Gene_expression.YRI, 
                             streamed_genotypes_matrix = binary_genotypes)


# Generate plots as above, but ONLY for associations with gene symbols
eQTL_to_barplots_scatterplot(output.dir = output.dir, 
                             output.filename = 'YRI_L1_eQTLs_Primary_Candidates_SNP1_Symbols_only_', 
                             plot_rows = 3, 
                             plot_columns = 3, 
                             indices_to_plot = c(1:1), 
                             combined_cis_trans_QTL_df = combined_cis_trans, 
                             all_snp_positions = all_snp_positions.YRI, 
                             gene_expression_df = Gene_expression.YRI, 
                             streamed_genotypes_matrix = binary_genotypes)


eQTL_to_barplots_scatterplot(output.dir = output.dir, 
                             output.filename = 'YRI_L1_eQTLs_Primary_Candidates_SNP2_Symbols_only_', 
                             plot_rows = 3, 
                             plot_columns = 3, 
                             indices_to_plot = c(4:4), 
                             combined_cis_trans_QTL_df = combined_cis_trans, 
                             all_snp_positions = all_snp_positions.YRI, 
                             gene_expression_df = Gene_expression.YRI, 
                             streamed_genotypes_matrix = binary_genotypes)
    

        
       
      
# PLOT GENOTYPE-GENE-TE RELATIONSHIPS * FOR TOP EUR TRIOS THAT DID NOT REPLICATE *   

# Update the snp details file by assigning RsIDs to the rownames
rownames(all_snp_positions.YRI) <- all_snp_positions.YRI$RsID    

# Make df to hold the relationships of interest
non.significant.trios <- combined_cis_trans[c(1:5), ]
non.significant.trios[,] <- NA

# Assign trio RsIDs
non.significant.trios[c(1:5), 'RsID'] <- c('rs11635336', 'rs11635336', 'rs9270493', 'rs9272300', 'rs9272222')

# Update trio snpid
non.significant.trios[c(1:5), 'snps'] <- all_snp_positions.YRI[non.significant.trios$RsID, 'snpid']

# Update trio gene
non.significant.trios[c(1:5), 'gene'] <- c('ENSG00000172345', 'ENSG00000172349', 'ENSG00000204308', 'ENSG00000237080', 'ENSG00000204315')

# Update trio TE
non.significant.trios[c(1:5), 'TE'] <- c('L1P4a:L1:LINE', 'L1P4a:L1:LINE', 'L1MEb:L1:LINE', 'L1MEb:L1:LINE', 'L1MEb:L1:LINE')

# Update trio symbols
non.significant.trios[c(1:5), 'symbol'] <- c('STARD5', 'IL16', 'RNF5', 'EHMT2-AS1', 'FKBPL')
 
# Update the snp details file by assigning my custom snpid to the rownames
rownames(all_snp_positions.YRI) <- all_snp_positions.YRI$snpid    

# Plot results
eQTL_to_barplots_scatterplot(output.dir = output.dir, 
                             output.filename = 'YRI_L1_eQTLs_NO_REPLICATION_Set1_', 
                             plot_rows = 3, 
                             plot_columns = 3, 
                             indices_to_plot = c(1:3), 
                             combined_cis_trans_QTL_df = non.significant.trios, 
                             all_snp_positions = all_snp_positions.YRI, 
                             gene_expression_df = Gene_expression.YRI, 
                             streamed_genotypes_matrix = binary_genotypes)      

# Plot results
eQTL_to_barplots_scatterplot(output.dir = output.dir, 
                             output.filename = 'YRI_L1_eQTLs_NO_REPLICATION_Set2_', 
                             plot_rows = 3, 
                             plot_columns = 3, 
                             indices_to_plot = c(4:5), 
                             combined_cis_trans_QTL_df = non.significant.trios, 
                             all_snp_positions = all_snp_positions.YRI, 
                             gene_expression_df = Gene_expression.YRI, 
                             streamed_genotypes_matrix = binary_genotypes)  
      
      
      
      
# Clean the environment
rm(list=ls())   
  
  
