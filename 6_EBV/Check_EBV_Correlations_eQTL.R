# Set strings as factors
options(stringsAsFactors = F)

# Load libraries  


# Load functions associated with this script.





# COMPARE EXPRESSION OF CANDIDATES WITH EBV


# Define output dir
output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/6_Integrate_cis_trans_eQTLs/EBV_Correlation_Plots/'

# Load Gene/TE expression (pre and post batch correction)
Gene_expression <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_eQTL/All_counts_EUR_358_filtered_VST.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')
Gene_expression.post <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_eQTL/All_counts_EUR_358_filtered_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Define expression values (pre correction)
exp.IL16 <- as.numeric(Gene_expression['ENSG00000172349', ])
exp.EBV <- as.numeric(Gene_expression['EBV_gene', ])

# Define expression values (POST correction)
exp.post.IL16 <- as.numeric(Gene_expression.post['ENSG00000172349', ])
exp.post.EBV <- as.numeric(Gene_expression.post['EBV_gene', ])

# Plot IL16 vs EBV (pre correction)
pdf(paste(output.dir, 'IL16_vs_EBV_Pre-batch',".pdf", sep=""), width = 8, height = 8)

plot(x = exp.IL16, 
     y = exp.EBV,
     xlab = 'Normalized IL16 log2(counts)',
     xlim = c(9, 13),
     ylab = 'Normalized EBV log2(counts)',
     ylim = c(12, 20),
     main = 'Pre-batch Correction Correlations',
     pch = 16,
     cex = 1)

dev.off()



# Plot IL16 vs EBV (post correction)
pdf(paste(output.dir, 'IL16_vs_EBV_post-batch',".pdf", sep=""), width = 8, height = 8)

plot(x = exp.post.IL16, 
     y = exp.post.EBV,
     xlab = 'Normalized IL16 log2(counts)',
     xlim = c(9, 13),
     ylab = 'Normalized EBV log2(counts)',
     ylim = c(12, 20),
     main = 'Post-batch Correction Correlations',
     pch = 16,
     cex = 1)

dev.off()

