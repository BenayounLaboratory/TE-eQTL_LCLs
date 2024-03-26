# Set strings as factors
options(stringsAsFactors = F)

# Load libraries  
library(CMplot) # Manhattan plots

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/6_Integrate_cis_trans_eQTLs/Integrate_eQTLs_and_Generate_Plots_Functions.R")     

# Define output directory
output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/6_Integrate_cis_trans_eQTLs/Manhattan_Plots/'





# MANHATTAN PLOTS



# Define input/general parameters

# Load snp_info
all_snp_positions.EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/RESOURCE_SNP_info_EUR.txt", header = T, stringsAsFactors = F, sep = '\t')

# Load table with empirical pval/FDR thresholds
empirical.thresholds <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/Significant_SNVs/EMPIRICAL_PVALUE_AND_FDR_THRESHOLDS.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Define FDR thresholds
BH.FDR.threshold <- 0.05

    

# Plot: EUR L1 trans eQTL

# Read external files. Only keep stats.
eQTL.stats <- readRDS('/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/Raw_eQTL_Results/L1_subfamily_trans_EUR.rds')
eQTL.stats <- eQTL.stats$all$eqtls

# Only plot results with pvalue smaller than 1e-3 (to reduce plotting time and filesize)
eQTL.stats <- eQTL.stats[which(eQTL.stats$pvalue <= 1e-3), ]

# Extract pvalue corresponding to the desired BH FDR
BH.pval_threshold <- max(eQTL.stats[eQTL.stats$FDR < BH.FDR.threshold, 'pvalue']) 
BH.pval_threshold # 2.31488e-08

# Define empirical FDR threshold pvalue
empirical.pval_threshold <- empirical.thresholds[c('EUR_L1_trans'), c('P_threshold')]

# Transform to CMplot format
CMplot_data <- eQTL_to_CMplot(eQTLs_for_plot = eQTL.stats, snp_pos = all_snp_positions.EUR)

# Define significant snps (to label on the Manhattan plot)
sig.snps <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/6_Integrate_cis_trans_eQTLs/Annotated_and_filtered_eQTLs/EUR_L1_eQTLs_annotated_1_filters_geneTEcorrelation_2023-05-22.txt", header = T, row.names = NULL, stringsAsFactors = F, sep = '\t')

    # For entries with no symbol, use the EnsemblID
    sig.snps[which(sig.snps$symbol == ""), 'symbol'] <- sig.snps[which(sig.snps$symbol == ""), 'gene']
    
    # Remove entries with duplicate EnsemblID
    sig.snps <- sig.snps[!duplicated(sig.snps$symbol), ]
    
    # Aggregate symbol column annotations by snpID
    aggregate.annotation <- aggregate(symbol ~ snps, sig.snps, FUN = toString)
    rownames(aggregate.annotation) <- aggregate.annotation$snps
    
    # To be able to obtain a CMPlot index, extract the significant results from the CMPlot object. NOTE: the threshold pvalues are included since their FDR is < 0.05
    CMPlot.filtered <- CMplot_data[which(CMplot_data$pvalue <= min(BH.pval_threshold,empirical.pval_threshold )), ]
    
    # Assign the 'snps' column as the rownames to the filtered CMPlot object (NOTE: this only works if the significant SNPs only appear once in the filtered list, which they do here)
    rownames(CMPlot.filtered) <- CMPlot.filtered$snps
    
    # To the annotated results, add column for the CMPlot index corresponding to each snp
    aggregate.annotation$CMPlot_index <- CMPlot.filtered[aggregate.annotation$snps, 'index']
    
# Make plot and save to file
pdf(file = paste(output.dir, "Plot_Manhattan_L1_trans-eQTLs_FDR5_EUR_", Sys.Date(), '.pdf', sep=""), width = 10, height = 5)
par(mar = c(5.1,6.1,4.1,2.1)) 

    CMplot(CMplot_data[, -c(5,6)],
           band = 0,
           plot.type = "m", 
           type = "p",
           col = c("black","grey60"),
           LOG10 = TRUE,
           amplify = FALSE,
           lwd.axis = 0.75,
           threshold.lwd = 0.75,
           threshold = c(BH.pval_threshold, empirical.pval_threshold),
           threshold.lty = c(1, 3),
           threshold.col = c('red', 'red'),
           ylim = c(3, 30),
           highlight = aggregate.annotation$CMPlot_index,
           highlight.col = NULL,
           highlight.cex = 0.25,
           highlight.text = aggregate.annotation$symbol,
           highlight.text.cex = 0.60,
           highlight.text.xadj = c(1, -1, -1, -1, 1, -1, -1),
           highlight.text.yadj = c(1, 0, 0, 1, 0, 1, 0),
           cex = 0.25, 
           cex.lab = 1,
           cex.axis = 1,
           chr.labels.angle = 45,
           #width = 20,
           #height = 10,
           #file = "tiff",
           verbose = TRUE,
           #dpi = 300,
           file.output = F)
dev.off()


# Loop over permutation files and plot

# Define output directory
perm.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/6_Integrate_cis_trans_eQTLs/Manhattan_Plots/L1_subfamily_trans_EUR_Permutations/'

# Define permutation data files
permutations.dir <- c("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/Raw_eQTL_Results/L1_subfamily_trans_EUR_Permutations/")
permutations.list <- list.files(permutations.dir, "\\.rds$", recursive = FALSE, full.names = TRUE)

# Define vector with permutation number (note that files are loaded alphabetically and not from smallest so largest number)
permutation_labels <- sort(as.character(1:length(permutations.list)))

# Define loop counter
counter <- 0

for (permutation_file in permutations.list) {
  
    # Update counter 
    counter <- counter + 1
    
    # Extract current permutation label
    ith_permutation_label <- permutation_labels[counter]
  
    # Read external files. Only keep stats.
    eQTL.stats <- readRDS(permutation_file)
    eQTL.stats <- eQTL.stats$all$eqtls
    
    # Only plot results with pvalue smaller than 1e-3 (to reduce plotting time and filesize)
    eQTL.stats <- eQTL.stats[which(eQTL.stats$pvalue <= 1e-3), ]
    
    # Transform to CMplot format
    CMplot_data <- eQTL_to_CMplot(eQTLs_for_plot = eQTL.stats, snp_pos = all_snp_positions.EUR)
    
    # Make plots and save to file
    pdf(file = paste(perm.output.dir, "Plot_Manhattan_L1_trans-eQTLs_FDR5_EUR_PERMUTATION_", ith_permutation_label, '_', Sys.Date(), '.pdf', sep=""), width = 10, height = 5)
    par(mar = c(5.1,6.1,4.1,2.1))
        CMplot(CMplot_data[, -c(5,6)],
               band = 0,
               plot.type = "m",
               type = "p",
               col = c("black","grey60"),
               LOG10 = TRUE,
               amplify = FALSE,
               lwd.axis = 0.75,
               threshold.lwd = 0.75,
               threshold = c(BH.pval_threshold, empirical.pval_threshold),
               threshold.lty = c(1, 3),
               threshold.col = c('red', 'red'),
               ylim = c(3, 30),
               #highlight = highlight_indices[1:20],
               #highlight.col = NULL,
               #highlight.cex = 0.25,
               #highlight.text = highlight_genes[1:20],
               #highlight.text.cex = 0.75,
               cex = 0.25,
               cex.lab = 1,
               cex.axis = 1,
               chr.labels.angle = 45,
               #width = 20,
               #height = 10,
               #file = "tiff",
               verbose = TRUE,
               #dpi = 300,
               file.output = F)
    dev.off()
}

        
           



        
   
    
        
# Plot: EUR Gene cis eQTL

# Read external files. Only keep stats.
eQTL.stats <- readRDS('/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/Raw_eQTL_Results/Gene_cis_EUR.rds')
eQTL.stats <- eQTL.stats$cis$eqtls

# Only plot results with pvalue smaller than 1e-3 (to reduce plotting time and filesize)
eQTL.stats <- eQTL.stats[which(eQTL.stats$pvalue <= 1e-3), ]

# Extract pvalue corresponding to user-defined FDR
BH.pval_threshold <- max(eQTL.stats[eQTL.stats$FDR < BH.FDR.threshold, 'pvalue']) 
BH.pval_threshold # 0.0004726647

# Define empirical FDR pvalue
empirical.pval_threshold <- empirical.thresholds[c('EUR_Gene_cis'), c('P_threshold')]

# Transform to CMplot format
CMplot_data <- eQTL_to_CMplot(eQTLs_for_plot = eQTL.stats, snp_pos = all_snp_positions.EUR)

# Make plots and save to file
pdf(file = paste(output.dir, "Plot_Manhattan_Gene_cis-eQTLs_FDR5_EUR_", Sys.Date(), '.pdf', sep=""), width = 10, height = 5)
par(mar = c(5.1,6.1,4.1,2.1)) 
    CMplot(CMplot_data[, -c(5,6)],
           plot.type = "m", 
           type = "p",
           col = c("black","grey60"),
           LOG10 = TRUE,
           amplify = FALSE,
           lwd.axis = 0.75,
           threshold.lwd = 0.75,
           threshold = c(BH.pval_threshold, empirical.pval_threshold),
           threshold.lty = c(1, 3),
           threshold.col = c('red', 'red'),
           ylim = c(3, 112),
           #highlight = highlight_indices[1:20],
           #highlight.col = NULL,
           #highlight.cex = 0.25,
           #highlight.text = highlight_genes[1:20],
           #highlight.text.cex = 0.75,
           cex = 0.25, 
           cex.lab = 1,
           cex.axis = 1,
           chr.labels.angle = 45,
           band = 0,
           #width = 20,
           #height = 10,
           #file = "tiff",
           verbose = TRUE,
           #dpi = 300,
           file.output = F)
dev.off()


# Loop over permutation files and plot

# Define output directory
perm.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/6_Integrate_cis_trans_eQTLs/Manhattan_Plots/Gene_cis_EUR_Permutations/'

# Define permutation data files
permutations.dir <- c("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/Raw_eQTL_Results/Gene_cis_EUR_Permutations/")
permutations.list <- list.files(permutations.dir, "\\.rds$", recursive=FALSE, full.names=TRUE)

# Define vector with permutation number (they're loaded alphabetically)
permutation_labels <- sort(as.character(1:length(permutations.list)))

# Define loop counter
counter <- 0

for (permutation_file in permutations.list) {
  
    # Update counter 
    counter <- counter + 1
    
    # Extract current permutation label
    ith_permutation_label <- permutation_labels[counter]
  
    # Read external files. Only keep stats.
    eQTL.stats <- readRDS(permutation_file)
    eQTL.stats <- eQTL.stats$cis$eqtls
    
    # Only plot results with pvalue smaller than 1e-3 (to reduce plotting time and filesize)
    eQTL.stats <- eQTL.stats[which(eQTL.stats$pvalue <= 1e-3), ]
    
    # Transform to CMplot format
    CMplot_data <- eQTL_to_CMplot(eQTLs_for_plot = eQTL.stats, snp_pos = all_snp_positions.EUR)
    
    # Make plots and save to file
    pdf(file = paste(perm.output.dir, "Plot_Manhattan_Gene_cis-eQTLs_FDR5_EUR_PERMUTATION_", ith_permutation_label, '_', Sys.Date(), '.pdf', sep=""), width = 10, height = 5)
    par(mar = c(5.1,6.1,4.1,2.1)) 
        CMplot(CMplot_data[, -c(5,6)],
               plot.type = "m", 
               type = "p",
               col = c("black","grey60"),
               LOG10 = TRUE,
               amplify = FALSE,
               lwd.axis = 0.75,
               threshold.lwd = 0.75,
               threshold = c(BH.pval_threshold, empirical.pval_threshold),
               threshold.lty = c(1, 3),
               threshold.col = c('red', 'red'),
               ylim = c(3, 112),
               #highlight = highlight_indices[1:20],
               #highlight.col = NULL,
               #highlight.cex = 0.25,
               #highlight.text = highlight_genes[1:20],
               #highlight.text.cex = 0.75,
               cex = 0.25, 
               cex.lab = 1,
               cex.axis = 1,
               chr.labels.angle = 45,
               band = 0,
               #width = 20,
               #height = 10,
               #file = "tiff",
               verbose = TRUE,
               #dpi = 300,
               file.output = F)
    dev.off()
  
}










# Plot: EUR L1 trans eQTL (taking into account PEER)

# Read external files. Only keep stats.
eQTL.stats <- readRDS('/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/Raw_eQTL_Results/L1_subfamily_trans_with_PEER_EUR.rds')
eQTL.stats <- eQTL.stats$all$eqtls

# Only plot results with pvalue smaller than 1e-3 (to reduce plotting time and filesize)
eQTL.stats <- eQTL.stats[which(eQTL.stats$pvalue <= 1e-3), ]

# Extract pvalue corresponding to the desired BH FDR
BH.pval_threshold <- max(eQTL.stats[eQTL.stats$FDR < BH.FDR.threshold, 'pvalue']) 
BH.pval_threshold # 2.444816e-07

# Transform to CMplot format
CMplot_data <- eQTL_to_CMplot(eQTLs_for_plot = eQTL.stats, snp_pos = all_snp_positions.EUR)
    
# Make plot and save to file
pdf(file = paste(output.dir, "Plot_Manhattan_L1_trans-eQTLs_with_PEER_FDR5_EUR_", Sys.Date(), '.pdf', sep=""), width = 10, height = 5)
par(mar = c(5.1,6.1,4.1,2.1)) 

    CMplot(CMplot_data[, -c(5,6)],
           band = 0,
           plot.type = "m", 
           type = "p",
           col = c("black","grey60"),
           LOG10 = TRUE,
           amplify = FALSE,
           lwd.axis = 0.75,
           threshold.lwd = 0.75,
           threshold = c(BH.pval_threshold),
           threshold.lty = c(1, 3),
           threshold.col = c('red', 'red'),
           ylim = c(3, 50),
           #highlight = aggregate.annotation$CMPlot_index,
           highlight.col = NULL,
           highlight.cex = 0.25,
           #highlight.text = aggregate.annotation$symbol,
           highlight.text.cex = 0.60,
           #highlight.text.xadj = c(1, -1, -1, -1, 1, -1, -1),
           #highlight.text.yadj = c(1, 0, 0, 1, 0, 1, 0),
           cex = 0.25, 
           cex.lab = 1,
           cex.axis = 1,
           chr.labels.angle = 45,
           #width = 20,
           #height = 10,
           #file = "tiff",
           verbose = TRUE,
           #dpi = 300,
           file.output = F)
dev.off()




        
      
    
    
# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/6_Integrate_cis_trans_eQTLs/Session_Info/'
    
    
sink(file = paste(dir.session_info,"Session_Info_Manhattan_Plots_with_PEER.txt", sep =""))
sessionInfo()
sink()      





# Clean the environment
rm(list=ls())       


