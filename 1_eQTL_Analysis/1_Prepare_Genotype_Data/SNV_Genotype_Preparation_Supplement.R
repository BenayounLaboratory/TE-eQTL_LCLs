## Generate snp position/mapping file and SNPID labels for the matrixeQTL 0/1/2 matrix
      
# Specify output dir
dir.genotypes.YRI <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_YRI_Genotypes/'
dir.genotypes.EUR <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/'

# Load initial data from plink BIM file
snp_info.YRI <- read.csv(paste('/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_YRI_Genotypes/Filtered_YRI_BOTH_SEXES_chr1-22.recode.plink.bim', sep = ''), header = F, stringsAsFactors = F, sep = '\t')
snp_info.EUR <- read.csv(paste('/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink.bim', sep = ''), header = F, stringsAsFactors = F, sep = '\t')

# Rename columns. !!!!!!!!!!!!!!!!!!! NOTE: If -keep-allele-order is used in plink, REF allele is automatically assigned to A2 (last column in BIM file). NOTE: The REF allele may or may not be the Major allele.
colnames(snp_info.YRI) <- c('chr', 'RsID', 'Morgan_Pos', 'pos', 'ALT', 'REF')
colnames(snp_info.EUR) <- c('chr', 'RsID', 'Morgan_Pos', 'pos', 'ALT', 'REF')

# Check if 'RsID' labels are unique
length(unique(snp_info.YRI$RsID)) == nrow(snp_info.YRI) # 1 mapping conflict. Conflict in dbSNP, not my workflow.
length(unique(snp_info.EUR$RsID)) == nrow(snp_info.EUR)

# Add snpid column. Format: row#_chr_pos. Important since some snps fall in the same position
snp_info.YRI$snpid <- paste(rownames(snp_info.YRI), snp_info.YRI$chr, snp_info.YRI$pos, sep = '_')
snp_info.EUR$snpid <- paste(rownames(snp_info.EUR), snp_info.EUR$chr, snp_info.EUR$pos, sep = '_')

# Rearrange columns, removing "Morgan position"
snp_info.YRI <- snp_info.YRI[, c('snpid', 'chr', 'pos', 'RsID', 'REF', 'ALT')]
snp_info.EUR <- snp_info.EUR[, c('snpid', 'chr', 'pos', 'RsID', 'REF', 'ALT')]

    # Save table
    write.table(snp_info.YRI, file = paste(dir.genotypes.YRI, "RESOURCE_SNP_info_YRI", ".txt", sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
    write.table(snp_info.EUR, file = paste(dir.genotypes.EUR, "RESOURCE_SNP_info_EUR", ".txt", sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)

# output the tranposed SNP ID column, which will be added to the genotype 012 matrix.
snpIDs_as_row_YRI <- t(snp_info.YRI[, 'snpid', drop = F])
snpIDs_as_row_EUR <- t(snp_info.EUR[, 'snpid', drop = F])

    # Save snpids
    write.table(snpIDs_as_row_YRI, file = paste(dir.genotypes.YRI, 'INPUT_Dependency_SNP_IDs_row_YRI', ".txt", sep =""), sep = " ", row.names = T, col.names = F, quote=F)
    write.table(snpIDs_as_row_EUR, file = paste(dir.genotypes.EUR, 'INPUT_Dependency_SNP_IDs_row_EUR', ".txt", sep =""), sep = " ", row.names = T, col.names = F, quote=F)