# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.


# Load libraries


# directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/5_Repeat_Family_Comparisons/Plots/'




# PREPARE DATA

# Load gene sets
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Gene_Set_Collections_for_GSEA.R')

    # Assign repeat subfamily to rownames of repeat geneset
    rownames(Repeat_family.gs) <- Repeat_family.gs$gene
    
    # Define repeat families used with GSEA
    GSEA.families <- c('L1 subfamilies', 'ERV1 subfamilies', 'ERVK subfamilies', 'ERVL subfamilies', 'ERVL-MaLR subfamilies', 'Alu subfamilies', 'hAT-Charlie subfamilies', 'TcMar-Tigger subfamilies', 'hAT-Tip100 subfamilies', 'Gypsy subfamilies')
    
    # Only keep repeat families used with GSEA
    Repeat_family.gs <- Repeat_family.gs[which(Repeat_family.gs$gs_name %in% GSEA.families), ]
    
# Make dataframe for each DESeq comparison
DESeq.Condition1 <- Repeat_family.gs
DESeq.Condition2 <- Repeat_family.gs

    # Add columns to hold log2 FC for each comparison
    DESeq.Condition1$Change <- NA
    DESeq.Condition2$Change <- NA

# Load IL16 OE (Condition 1) DESeq results
DESeq.res <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/3_DESeq/DESeq_Results_GM12878_OE/All_Genes_IL16.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

    # Obtain list of expressed repeats
    my.repeats <- Repeat_family.gs[Repeat_family.gs$gene %in% rownames(DESeq.res), 'gene']

    # Add the log2FC to the df made earlier
    DESeq.Condition1 <- DESeq.Condition1[my.repeats, ]
    DESeq.Condition1$Change <- DESeq.res[my.repeats, 'log2FoldChange']
    
# Load STARD5 OE (Condition 2) DESeq results
DESeq.res <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/3_DESeq/DESeq_Results_GM12878_OE/All_Genes_STARD5.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

    # Obtain list of expressed repeats
    my.repeats <- Repeat_family.gs[Repeat_family.gs$gene %in% rownames(DESeq.res), 'gene']

    # Add the log2FC to the df made earlier
    DESeq.Condition2 <- DESeq.Condition2[my.repeats, ]
    DESeq.Condition2$Change <- DESeq.res[my.repeats, 'log2FoldChange']
    

# Organize family names alphabetically
DESeq.Condition1 <- DESeq.Condition1[order(DESeq.Condition1$gs_name), ]
DESeq.Condition2 <- DESeq.Condition2[order(DESeq.Condition2$gs_name), ]

# Lock in factors
DESeq.Condition1$gs_name <- factor(DESeq.Condition1$gs_name , levels = GSEA.families)
DESeq.Condition2$gs_name <- factor(DESeq.Condition2$gs_name , levels = GSEA.families)

    
    
    


    

# RUN STATISTICS (ONE SAMPLE WILCOXON TEST)

# Make a function to calculate pvalues
calc_pval <- function(input.DESeq.Condition){
  
      # Make a dataframe to hold pvalues
      all.my.stats <- data.frame(Families = GSEA.families, pval = NA, FDR = NA)
      rownames(all.my.stats) <- GSEA.families
      
      # Loop over each families and calculate pvalues
      for (TE_family in all.my.stats$Families) {
        
          # Extract the log2FC for the ith repeat family
          ith_changes <- input.DESeq.Condition[which(input.DESeq.Condition$gs_name == TE_family), 'Change']
        
          # Run Wilcoxon on ith log2FC
          ith_stats <- wilcox.test(ith_changes, mu = 0, alternative = "two.sided", exact = TRUE)
          
          # Fill in the pvalue in the df
          all.my.stats[TE_family, 'pval'] <- ith_stats$p.value
          
      }
      
      # Apply FDR correction to pvalues
      all.my.stats$FDR <- p.adjust(all.my.stats$pval, method = 'fdr')
      
      # Only keep 3 sig figures for FDR
      all.my.stats$FDR_short <- signif(all.my.stats$FDR, 3)
      
      # Return the results
      return(all.my.stats)
  
} # END FUNCTION


# Calculate pvalues for each condition
stats.condition1 <- calc_pval(input.DESeq.Condition = DESeq.Condition1)
stats.condition2 <- calc_pval(input.DESeq.Condition = DESeq.Condition2)








# GENERATE PLOTS
my.label.yval <- -0.85

# IL16 OE
pdf(paste(my.output.dir, "Boxplot_Repeat_log2FC_OE_IL16", ".pdf", sep=""), width = 3.25, height = 4)

    # Remove white space from header (3rd value in vector)
    par(mar=c(5,4,1,2)+0.1)

    # make boxplots for genotype segregated gene expression
    boxplot(Change ~ gs_name, 
            data = DESeq.Condition1,
            ylim = c(-1.0, 1.0),
            outline = F, 
            #col = boxplot.col,
            main = c('IL16 OE Effects'),
            cex.main = 0.70,
            ylab = c('log2(Fold Change)'),
            xlab = c(''),
            cex.lab = 0.7,
            cex.axis = 0.7,
            xaxt = 'n')
    
    # overlay the points
    stripchart(Change ~ gs_name, data = DESeq.Condition1, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 0.3, col = 'red') 
    
    # Add line at y = 0
    abline(h = 0, lty = 5, col = 'black')
    
    # Specify axis labels
    axis(1, at = 1:10, labels = GSEA.families, las = 2, cex.axis = 0.70)

    # Add pvalues
    text(1, my.label.yval, stats.condition1[1, 'FDR_short'], cex = 0.7, srt = 90)
    text(2, my.label.yval, stats.condition1[2, 'FDR_short'], cex = 0.7, srt = 90)
    text(3, my.label.yval, stats.condition1[3, 'FDR_short'], cex = 0.7, srt = 90)
    text(4, my.label.yval, stats.condition1[4, 'FDR_short'], cex = 0.7, srt = 90)
    text(5, my.label.yval, stats.condition1[5, 'FDR_short'], cex = 0.7, srt = 90)
    text(6, my.label.yval, stats.condition1[6, 'FDR_short'], cex = 0.7, srt = 90)
    text(7, my.label.yval, stats.condition1[7, 'FDR_short'], cex = 0.7, srt = 90)
    text(8, my.label.yval, stats.condition1[8, 'FDR_short'], cex = 0.7, srt = 90)
    text(9, my.label.yval, stats.condition1[9, 'FDR_short'], cex = 0.7, srt = 90)
    text(10, my.label.yval, stats.condition1[10, 'FDR_short'], cex = 0.7, srt = 90)

    
# End pdf
dev.off()    
    


# STARD5 OE
my.label.yval <- -0.85

pdf(paste(my.output.dir, "Boxplot_Repeat_log2FC_OE_STARD5", ".pdf", sep=""), width = 3.25, height = 4)

    # Remove white space from header (3rd value in vector)
    par(mar=c(5,4,1,2)+0.1)

    # make boxplots for genotype segregated gene expression
    boxplot(Change ~ gs_name, 
            data = DESeq.Condition2,
            ylim = c(-1.0, 1.0),
            outline = F, 
            #col = boxplot.col,
            main = c('STARD5 OE Effects'),
            cex.main = 0.70,
            ylab = c('log2(Fold Change)'),
            xlab = c(''),
            cex.lab = 0.7,
            cex.axis = 0.7,
            xaxt = 'n')
    
    # overlay the points
    stripchart(Change ~ gs_name, data = DESeq.Condition2, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 0.3, col = 'red') 
    
    # Add line at y = 0
    abline(h = 0, lty = 5, col = 'black')
    
    # Specify axis labels
    axis(1, at = 1:10, labels = GSEA.families, las = 2, cex.axis = 0.70)

    # Add pvalues
    text(1, my.label.yval, stats.condition2[1, 'FDR_short'], cex = 0.7, srt = 90)
    text(2, my.label.yval, stats.condition2[2, 'FDR_short'], cex = 0.7, srt = 90)
    text(3, my.label.yval, stats.condition2[3, 'FDR_short'], cex = 0.7, srt = 90)
    text(4, my.label.yval, stats.condition2[4, 'FDR_short'], cex = 0.7, srt = 90)
    text(5, my.label.yval, stats.condition2[5, 'FDR_short'], cex = 0.7, srt = 90)
    text(6, my.label.yval, stats.condition2[6, 'FDR_short'], cex = 0.7, srt = 90)
    text(7, my.label.yval, stats.condition2[7, 'FDR_short'], cex = 0.7, srt = 90)
    text(8, my.label.yval, stats.condition2[8, 'FDR_short'], cex = 0.7, srt = 90)
    text(9, my.label.yval, stats.condition2[9, 'FDR_short'], cex = 0.7, srt = 90)
    text(10, my.label.yval, stats.condition2[10, 'FDR_short'], cex = 0.7, srt = 90)

    
# End pdf
dev.off()        
    
    
    
    
    




    
# Clean the environment
rm(list=ls())
    
    
