# Set strings as factors
options(stringsAsFactors = F)

# Define output directory
dir.output <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/Structural_Variant_Counts/'

# Load 0/1/2 matrices for L1/Alu insertions
L1_all_Insertions <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/Structural_Variant_Counts/LINE1_Insertions/ALL_LINE1_insertions_transposed.recode.plink.012", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')
Alu_all_Insertions <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/Structural_Variant_Counts/ALU_Insertions/ALL_ALU_insertions_transposed.recode.plink.012", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Load 0/1/2 matrices for L1/Alu deletions
L1_all_deletions <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/Structural_Variant_Counts/LINE1_deletions/ALL_LINE1_deletions_transposed.recode.plink.012", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')
Alu_all_deletions <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/Structural_Variant_Counts/ALU_deletions/ALL_ALU_deletions_transposed.recode.plink.012", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Make table to hold summary insertion/deletion info
TE_net_inserts <- as.data.frame(matrix(0, nrow = 2504, ncol = 6))
colnames(TE_net_inserts) <- c('L1_inserts', 'L1_deletions', 'L1_net', 'Alu_inserts', 'Alu_deletions', 'Alu_net')
rownames(TE_net_inserts) <- colnames(L1_all_Insertions)

# Fill in: The total number of insertions detected
TE_net_inserts[colnames(L1_all_Insertions), 'L1_inserts'] <- colSums(L1_all_Insertions)
TE_net_inserts[colnames(Alu_all_Insertions), 'Alu_inserts'] <-colSums(Alu_all_Insertions)

# Fill in: The total number of deletions detected
TE_net_inserts[colnames(L1_all_deletions), 'L1_deletions'] <- colSums(L1_all_deletions)
TE_net_inserts[colnames(Alu_all_deletions), 'Alu_deletions'] <- colSums(Alu_all_deletions)

# Fill in: The net number of insertions relative to the reference (insertions - deletions)
TE_net_inserts$L1_net <- TE_net_inserts$L1_inserts - TE_net_inserts$L1_deletions
TE_net_inserts$Alu_net <- TE_net_inserts$Alu_inserts - TE_net_inserts$Alu_deletions

# Fill in: The net L1 + Alu insertions
TE_net_inserts$Net_Both <-  TE_net_inserts$L1_net + TE_net_inserts$Alu_net

# Save the table; net insertions may be used as a covariate in some analyses
write.table(TE_net_inserts, file = paste(dir.output, "L1_Alu_Insertions_Deletions_Per_Sample", '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)

# Remove unneeded variables
rm(list=ls())


  

