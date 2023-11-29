# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(eQTLMAPT) # To run mediation analysis
library(BEDMatrix) # To stream genotypes 

# Define output dir
output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/6_Mediation_Analysis/Results/'



# MEDIATION ANALYSIS WITH EQTLMAPT --------------------------



# Load SNVs that pass the tri-part analysis
combined_cis_trans.EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/4_Integrate_cis_trans_eQTLs/Annotated_and_filtered_eQTLs/EUR_L1_intergenic_distal_eQTLs_annotated_1_filters_geneTEcorrelation.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')
#combined_cis_trans.YRI <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/7_Replication_of_eQTLs_in_YRI/Annotated_and_filtered_eQTLs/YRI_L1_eQTLs_annotated_1_filters_geneTEcorrelation_2023-05-23.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')

    # Define unique SNVs
    unique.snps.EUR <- unique(combined_cis_trans.EUR$RsID)
    #unique.snps.YRI <- unique(combined_cis_trans.YRI$RsID)

# Define genotype inputs

    # Specify the genotype BED file
    plink_file_path.EUR <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink.bed'
    #plink_file_path.YRI <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_YRI_Genotypes/Filtered_YRI_BOTH_SEXES_chr1-22.recode.plink.bed'
    
    # Make R object to stream snp genotypes into memory (rows are samples and columns are snps). # NOTE: If used, MT snps will need to be assigned a snpid, since they dont have rsid.
    binary_genotypes.EUR <- BEDMatrix(path = plink_file_path.EUR, simple_names = T)
    #binary_genotypes.YRI <- BEDMatrix(path = plink_file_path.YRI, simple_names = T)

    # Make a vector of the sample names (to make sure all inputs are in the same order)
    sample.ids.EUR <- rownames(binary_genotypes.EUR)
    #sample.ids.YRI <- rownames(binary_genotypes.YRI)

    # Extract genotypes for snps to be tested
    snp.dat.EUR <- binary_genotypes.EUR[sample.ids.EUR, unique.snps.EUR]
    #snp.dat.YRI <- binary_genotypes.YRI[sample.ids.YRI, unique.snps.YRI]
    
    # Transpose genotypes to meet software input format
    snp.dat.EUR <- t(snp.dat.EUR)
    #snp.dat.YRI <- t(snp.dat.YRI)
    
    # Make snp.dat as data frame
    snp.dat.EUR <- as.data.frame(snp.dat.EUR)
    #snp.dat.YRI <- as.data.frame(snp.dat.YRI)
    
    # Make a dataframe to hold the index of each unique snp
    snp.indices.EUR <- data.frame(row.names = rownames(snp.dat.EUR),
                                  index = 1:nrow(snp.dat.EUR))
    
    # snp.indices.YRI <- data.frame(row.names = rownames(snp.dat.YRI),
    #                               index = 1:nrow(snp.dat.YRI))
    
# Load Gene/TE expression (use INT data with covariates removed)
INT.expression.EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/2_Prepare_RNASeq/Processed_counts_eQTL/TElocal_counts_EUR_358_filtered_aggregate_split_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV_INT_All_Genes_TEs.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')
#INT.expression.YRI <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_eQTL/All_counts_YRI_86_filtered_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV_INT_All_Genes_TEs.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

    # Reorder samples for uniformity across inputs
    INT.expression.EUR <- INT.expression.EUR[, sample.ids.EUR]
    #INT.expression.YRI <- INT.expression.YRI[, sample.ids.YRI]
    
    # Make a dataframe to hold the index of each gene/TE
    expression.indices.EUR <- data.frame(row.names = rownames(INT.expression.EUR),
                                         index = 1:nrow(INT.expression.EUR))
    
    # expression.indices.YRI <- data.frame(row.names = rownames(INT.expression.YRI),
    #                                      index = 1:nrow(INT.expression.YRI))
   
# Prepare the trio matrix

    # Extract the necessary components from my sig results table
    trios.EUR <- combined_cis_trans.EUR[, c('RsID', 'gene', 'TE')]
    #trios.YRI <- combined_cis_trans.YRI[, c('RsID', 'gene', 'TE')]
    
    # Update the RsID column with snp indices
    trios.EUR$RsID <- snp.indices.EUR[trios.EUR$RsID, 'index']
    #trios.YRI$RsID <- snp.indices.YRI[trios.YRI$RsID, 'index']
    
    # Update the gene column with gene indices
    trios.EUR$gene <- expression.indices.EUR[trios.EUR$gene, 'index']
    #trios.YRI$gene <- expression.indices.YRI[trios.YRI$gene, 'index']
    
    # Update the TE column with TE indices
    trios.EUR$TE <- expression.indices.EUR[trios.EUR$TE, 'index']
    #trios.YRI$TE <- expression.indices.YRI[trios.YRI$TE, 'index']

# Make a known confounder matrix 
    
    # If known confounders have been removed from the expression table being used, then assign a value of 1
    conf.EUR <- as.data.frame(matrix(1, nrow = 1, ncol = length(sample.ids.EUR)))
    colnames(conf.EUR) <- sample.ids.EUR
    
    #conf.YRI <- as.data.frame(matrix(1, nrow = 1, ncol = length(sample.ids.YRI)))
    #colnames(conf.YRI) <- sample.ids.YRI
  
# Run eqtlmapt in basic mode 

    # Generate a cluster with 7 nodes for parallel computing  
    my.cluster = makeCluster(7)  
    
    # Genomic mediation, adaptive permutations
    eqtlmapt.output.EUR <- gmap(snp.dat = as.matrix(snp.dat.EUR),
                                fea.dat = as.matrix(INT.expression.EUR),
                                conf = conf.EUR,
                                trios.idx = as.matrix(trios.EUR),
                                cl = my.cluster,
                                Minperm = 100,
                                Maxperm = 30000)
    
    # eqtlmapt.output.YRI <- gmap(snp.dat = as.matrix(snp.dat.YRI),
    #                             fea.dat = as.matrix(INT.expression.YRI),
    #                             conf = conf.YRI,
    #                             trios.idx = as.matrix(trios.YRI),
    #                             cl = my.cluster,
    #                             Minperm = 100,
    #                             Maxperm = 30000)
    
    # Stop the parallel cluster
    stopCluster(my.cluster)
    
    # Convert the results to a dataframe
    mediation.results.EUR <- as.data.frame(eqtlmapt.output.EUR)
    #mediation.results.YRI <- as.data.frame(eqtlmapt.output.YRI)
    
    # Add snp-gene-te info to mediation results
    mediation.results.EUR <- cbind(combined_cis_trans.EUR[, c(1:4)],
                                   mediation.results.EUR)
    
    # mediation.results.YRI <- cbind(combined_cis_trans.YRI[, c(1:4)],
    #                                mediation.results.YRI)
    
    # Use empirical pvalues to calculate BH FDR
    mediation.results.EUR$EMP_P_FDR <- p.adjust(mediation.results.EUR$empirical.p, method = 'fdr')
    #mediation.results.YRI$EMP_P_FDR <- p.adjust(mediation.results.YRI$empirical.p, method = 'fdr')
      
    # Save data
    write.table(mediation.results.EUR, file = paste(output.dir, "Mediation_Analysis_EUR_L1_intergenic_distal_Candidates_gmap_with_eQTL_data", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
    #write.table(mediation.results.YRI, file = paste(output.dir, "Mediation_Analysis_YRI_L1_Candidates_gmap_with_eQTL_data_", Sys.Date(), '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
    
    # Filter mediation results on empirical pvalue FDR
    mediation.results.EUR.sig <- mediation.results.EUR[which(mediation.results.EUR$EMP_P_FDR < 0.05), ]
    #mediation.results.YRI.sig <- mediation.results.YRI[which(mediation.results.YRI$EMP_P_FDR < 0.05), ]
    
        # Save data
        write.table(mediation.results.EUR.sig, file = paste(output.dir, "Mediation_Analysis_EUR_L1_intergenic_distal_Candidates_gmap_with_eQTL_data_FDR5", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
        # write.table(mediation.results.YRI.sig, file = paste(output.dir, "Mediation_Analysis_YRI_L1_Candidates_gmap_with_eQTL_data_FDR5_", Sys.Date(), '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F) # no results pass filtering
        
      
    

        

# PLOT PROPORTION MEDIATED  --------------------------



# Load mediation results
mediation.EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/6_Mediation_Analysis/Results/Mediation_Analysis_EUR_L1_intergenic_distal_Candidates_gmap_with_eQTL_data_FDR5.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')
#mediation.YRI <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/8_Mediation_Analysis/Results/Mediation_Analysis_YRI_L1_Candidates_gmap_with_eQTL_data_2023-05-07.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')

# Order results by nominal pval
mediation.EUR <- mediation.EUR[order(mediation.EUR$nominal.p, decreasing = FALSE),]
#mediation.YRI <- mediation.YRI[order(mediation.YRI$nominal.p, decreasing = FALSE),]
    
# Remove entries without a protein coding gene
mediation.EUR <- mediation.EUR[which(mediation.EUR$symbol != ""),]
#mediation.YRI <- mediation.YRI[which(mediation.YRI$symbol != ""),]
    
# Keep first instance of gene symbol
mediation.EUR <- mediation.EUR[!duplicated(mediation.EUR$symbol), ]   
#mediation.YRI <- mediation.YRI[!duplicated(mediation.YRI$symbol), ]   

# Order results by proportion mediated
mediation.EUR <- mediation.EUR[order(mediation.EUR$beta.change, decreasing = TRUE),]
#mediation.YRI <- mediation.YRI[order(mediation.YRI$beta.change, decreasing = TRUE),]

# Update FDR values to scientific notation
mediation.EUR$EMP_P_FDR <- formatC(mediation.EUR$EMP_P_FDR, format = "e", digits = 2)
#mediation.YRI$EMP_P_FDR <- formatC(mediation.YRI$EMP_P_FDR, format = "e", digits = 2)

# Define group levels (needed to specify order in the plot)
mediation.EUR$symbol <- factor(mediation.EUR$symbol, levels = mediation.EUR$symbol)
#mediation.YRI$symbol <- factor(mediation.YRI$symbol, levels = mediation.YRI$symbol)

# Generate EUR bar plots
pdf(paste(output.dir, "Plot_Most_Significant_Proportion_Mediated_EUR_intergenic_distal", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(mediation.EUR$beta.change ~ mediation.EUR$symbol,
            ylab = "Proportion of SNV-L1 Effect Mediated",
            ylim = c(-0.2, 0.6),
            xlab = "",
            las = 2,
            cex.names = 1
            )

    # Add FDR to plot
    text(0.7, 0.5, mediation.EUR[1, 'EMP_P_FDR'], cex = 1)
    
dev.off()


# # Generate YRI bar plots
# pdf(paste(output.dir, "Plot_Most_Significant_Proportion_Mediated_YRI_", Sys.Date(), ".pdf", sep=""), width = 6, height = 6)
# par(mar=c(6,4,4,1)+.1)
# 
#     # Plot proportion mediated as a barplot
#     barplot(mediation.YRI$beta.change ~ mediation.YRI$symbol,
#             ylab = "Proportion of SNV-L1 Effect Mediated",
#             ylim = c(0, 0.6),
#             xlab = "",
#             las = 2,
#             cex.names = 1
#             )
# 
#     # Add FDR to plot
#     text(0.7, 0.5, mediation.YRI[1, 'EMP_P_FDR'], cex = 1)
#     text(1.9, 0.5, mediation.YRI[2, 'EMP_P_FDR'], cex = 1)
# 
#     
# dev.off()





# PLOT DIRECT/INDIRECT/TOTAL EFFECTS  --------------------------



# Load mediation results
mediation.EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1b_eQTL_Analysis_TElocal/6_Mediation_Analysis/Results/Mediation_Analysis_EUR_L1_intergenic_distal_Candidates_gmap_with_eQTL_data_FDR5.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')

# Calculate Mediation/Indirect Effects (this is equal to beta.total - beta)
mediation.EUR$Indirect <- mediation.EUR$beta.total - mediation.EUR$beta

# Order results by nominal pval
mediation.EUR <- mediation.EUR[order(mediation.EUR$nominal.p, decreasing = FALSE),]

# Remove entries without a protein coding gene
mediation.EUR <- mediation.EUR[which(mediation.EUR$symbol != ""),]

# Update row indices
rownames(mediation.EUR) <- 1:nrow(mediation.EUR)

# Specify row numbers to plot (corresponding to SNV-gene interactions in the paper figures)
rows_to_plot.EUR <- c(1) 

# Extract data to plot
data_to_plot.EUR <- mediation.EUR[rows_to_plot.EUR, ]

# Assign gene names to rownames
rownames(data_to_plot.EUR) <- data_to_plot.EUR$symbol

# Remove columns we don't need
data_to_plot.EUR <- data_to_plot.EUR[, c('beta', 'Indirect', 'beta.total')]

# Update colnames to match plot labels
colnames(data_to_plot.EUR) <- c('Direct', 'Indirect', 'Total')

# Transpose the data to plot
data_to_plot.EUR <- as.data.frame(t(data_to_plot.EUR))

# Add a column with the beta label
data_to_plot.EUR$all_betas <- rownames(data_to_plot.EUR)

# Separately define genes with positive and negative betas
data_to_plot.EUR.pos <- data_to_plot.EUR[, c(1,2), drop = FALSE]
#data_to_plot.EUR.neg <- data_to_plot.EUR[, c(1,2), drop = FALSE]

# Define group levels (needed to specify order in the plot)
data_to_plot.EUR.pos$all_betas <- factor(data_to_plot.EUR.pos$all_betas, levels = data_to_plot.EUR.pos$all_betas)
#data_to_plot.EUR.neg$all_betas <- factor(data_to_plot.EUR.neg$all_betas, levels = data_to_plot.EUR.neg$all_betas)




pdf(paste(output.dir, "Barplot_Direct_Indirect_Total_Effects_EUR_Positive_Intergenic_Distal",".pdf", sep=""), width = 7, height = 7)
par(mfrow=c(3,3)) # Added to produce images of similar shape across analyses

    # Loop over genes
    for (ith_gene in 1:(ncol(data_to_plot.EUR.pos)-1)) {

      # Subset the ith gene column and the all_betas column
      ith_data <- data_to_plot.EUR.pos[, c(ith_gene, ncol(data_to_plot.EUR.pos))]

      # Change the gene symbol to a generic label
      colnames(ith_data)[1] <- 'Effect'

      # Plot proportion mediated as a barplot
      barplot(ith_data$Effect ~ ith_data$all_betas,
              ylab = "Effect Coefficient",
              ylim = c(0, 1),
              xlab = "",
              las = 2,
              cex.names = 1,
              )

    }


dev.off()


# pdf(paste(output.dir, "Barplot_Direct_Indirect_Total_Effects_EUR_Negative_Intergenic_Distal",".pdf", sep=""), width = 7, height = 7)
# par(mfrow=c(3,3)) # Added to produce images of similar shape across analyses
# 
#     # Loop over genes
#     for (ith_gene in 1:(ncol(data_to_plot.EUR.neg)-1)) {
#       
#       # Subset the ith gene column and the all_betas column
#       ith_data <- data_to_plot.EUR.neg[, c(ith_gene, ncol(data_to_plot.EUR.neg))]
#       
#       # Change the gene symbol to a generic label
#       colnames(ith_data)[1] <- 'Effect'
#       
#       # Plot proportion mediated as a barplot
#       barplot(ith_data$Effect ~ ith_data$all_betas,
#               ylab = "Effect Coefficient",
#               ylim = c(0, -1),
#               xlab = "",
#               las = 2,
#               cex.names = 1,
#               )
#       
#     }
# 
# 
# dev.off()






# Clean the environment
rm(list=ls())   
