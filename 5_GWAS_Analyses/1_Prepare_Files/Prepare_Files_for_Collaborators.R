# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(biomaRt) # for mapping RsIds to hg19 positions
    
# IMPORTANT NOTE: THE SCRIPT BELOW ASSUMES THAT THERE IS ONLY ONE L1 EQTL PVALUE/FDR FOR EACH SNV (WHICH THERE IS, IN MY ANALYSIS)
    
    

        
# Specify and output directory
output.dir <- c('/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/5_GWAS_Analyses/1_Prepare_Files/')    
    
# Load eQTL SNP lists
EUR.L1.snps <- read.csv(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/Significant_SNVs/L1_trans_eQTLs_SIGNIFICANT_EUR_2023-05-22.txt', header = TRUE, sep = "")
YRI.L1.snps <- read.csv(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/7_Replication_of_eQTLs_in_YRI/Scan_for_eQTLs/L1_subfamily_trans_eQTLs_YRI_FDR5_2023-05-23.txt', header = TRUE, sep = "")
    
    # Remove duplicate RsIds
    EUR.L1.snps <- EUR.L1.snps[!duplicated(EUR.L1.snps$RsID), ]
    YRI.L1.snps <- YRI.L1.snps[!duplicated(YRI.L1.snps$RsID), ]

# Load table with details of each SNP used in the eQTL analysis
my.SNP.metadata <- read.csv(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/RESOURCE_SNP_info_EUR.txt', header = TRUE, sep = "")

    # Only keep snps significant in the eQTL analysis (using my custom IDs, which deal with duplicate RsIDs)
    my.SNP.metadata <- my.SNP.metadata[which(my.SNP.metadata$snpid %in% EUR.L1.snps$snps), ]
    
    # Assign RsID as the rowname
    rownames(my.SNP.metadata) <- my.SNP.metadata$RsID

# Add REF/ALT allele info to the eQTL results table (note: I'm going to use the EUR table as a template since the YRI analyses are a subset of the EUR analyses)
EUR.L1.snps[, 'REF'] <- my.SNP.metadata[EUR.L1.snps$RsID, 'REF']
EUR.L1.snps[, 'ALT'] <- my.SNP.metadata[EUR.L1.snps$RsID, 'ALT']

# Add CHR and BP info to the eQTL results table
EUR.L1.snps[, 'HG38_CHR'] <- my.SNP.metadata[EUR.L1.snps$RsID, 'chr']
EUR.L1.snps[, 'HG38_BP'] <- my.SNP.metadata[EUR.L1.snps$RsID, 'pos']
    
# Make a new df to hold the final results that will be shared with collaborators
snp_info_for_collabs <- EUR.L1.snps

    # Assign RsID to rowname for easier manipulation
    rownames(snp_info_for_collabs) <- snp_info_for_collabs$RsID

    # Keep and reorganize only necessary columns
    snp_info_for_collabs <- snp_info_for_collabs[, c('RsID', 'HG38_CHR', 'HG38_BP', 'REF', 'ALT', 'pvalue', 'FDR')]
    
    # Update colnames for clarity
    colnames(snp_info_for_collabs) <- c('SNP', 'HG38_CHR', 'HG38_BP', 'REF', 'ALT', 'EUR_pval', 'EUR_FDR')
    
    # Add NA columns to fill in AFR pval/FDR if available
    snp_info_for_collabs[, 'AFR_pval'] = NA
    snp_info_for_collabs[, 'AFR_FDR'] = NA

# Fill in AFR pval and FDR
snp_info_for_collabs[YRI.L1.snps$RsID, 'AFR_pval'] <- YRI.L1.snps$pvalue
snp_info_for_collabs[YRI.L1.snps$RsID, 'AFR_FDR'] <- YRI.L1.snps$FDR

# Obtain hg19 coordinates for eQTL RsIds (NOTE: SNVs without an RsID will not map)
    
    # Select mart
    # ensembl.snp.mart <- useEnsembl("snp", dataset = "hsapiens_snp", GRCh = "37", version = 107, verbose = TRUE)
    ensembl.snp.mart <- useEnsembl("snp", dataset = "hsapiens_snp", GRCh = "37", verbose = TRUE)
    
    #get genomic position
    snp_hg19_mapping_table <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end"), 
                                    filters = "snp_filter", 
                                    values = snp_info_for_collabs$SNP, 
                                    mart = ensembl.snp.mart, 
                                    uniqueRows = TRUE,
                                    verbose = FALSE)

    # order the mapping table by chromosome number (this is done to move patch positions to the end of the list, and preferentially remove them if a non-patched position exists)
    snp_hg19_mapping_table <- snp_hg19_mapping_table[order(snp_hg19_mapping_table$chr_name, decreasing = FALSE), ]
    
    # Check names of duplicate RsIds
    duplicate_rsid_chromosomes <- unique(snp_hg19_mapping_table[duplicated(snp_hg19_mapping_table$refsnp_id), 'chr_name'])
    
    # remove duplicate RsId entries (should correspond to genome patches based on the previous line of code)
    snp_hg19_mapping_table <- snp_hg19_mapping_table[!duplicated(snp_hg19_mapping_table$refsnp_id), ]
    
# Add hg19 coordinates to GTAS table
    
    # Add the hg19 chromosome info to the table
    snp_info_for_collabs[snp_hg19_mapping_table$refsnp_id, 'HG19_CHR'] <- snp_hg19_mapping_table$chr_name
    
    # Add the hg19 base position info to the table
    snp_info_for_collabs[snp_hg19_mapping_table$refsnp_id, 'HG19_BP'] <- snp_hg19_mapping_table$chrom_start

# Reorganize the table
snp_info_for_collabs <- snp_info_for_collabs[, c(1:3, 10:11, 4:9)]

# Save the table
write.table(snp_info_for_collabs, file = paste(output.dir, "Bravo_eQTL_SNVs_for_Collaborators_", Sys.Date(), '.txt', sep=""), row.names = FALSE, col.names = TRUE, sep = "\t", na = "NA", quote = F)

    
        


# Clear variables from the environment
rm(list=ls())
    
    
        
