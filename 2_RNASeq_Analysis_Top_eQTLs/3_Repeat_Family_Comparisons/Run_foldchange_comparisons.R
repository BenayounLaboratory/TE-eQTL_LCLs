# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.


# Load libraries


# directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/3_Repeat_Family_Comparisons/Plots/'




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
DESeq.Condition3 <- Repeat_family.gs
DESeq.Condition4 <- Repeat_family.gs
DESeq.Condition5 <- Repeat_family.gs
    
    # Add columns to hold log2 FC for each comparison
    DESeq.Condition1$Change <- NA
    DESeq.Condition2$Change <- NA
    DESeq.Condition3$Change <- NA
    DESeq.Condition4$Change <- NA
    DESeq.Condition5$Change <- NA

# Load IL16 SNV (Condition 1) DESeq results
DESeq.res <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/1_DESeq_Genotypes/DESeq_Results/EUR_rs11635336_All_Genes.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

    # Obtain list of expressed repeats
    my.repeats <- Repeat_family.gs[Repeat_family.gs$gene %in% rownames(DESeq.res), 'gene']

    # Add the log2FC to the df made earlier
    DESeq.Condition1 <- DESeq.Condition1[my.repeats, ]
    DESeq.Condition1$Change <- DESeq.res[my.repeats, 'log2FoldChange']
    
# Load HLA SNV (Condition 2) DESeq results
DESeq.res <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/1_DESeq_Genotypes/DESeq_Results/EUR_rs9271894_All_Genes.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

    # Obtain list of expressed repeats
    my.repeats <- Repeat_family.gs[Repeat_family.gs$gene %in% rownames(DESeq.res), 'gene']

    # Add the log2FC to the df made earlier
    DESeq.Condition2 <- DESeq.Condition2[my.repeats, ]
    DESeq.Condition2$Change <- DESeq.res[my.repeats, 'log2FoldChange']
    
# Load HSD17B12 SNV (Condition 3) DESeq results
DESeq.res <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/1_DESeq_Genotypes/DESeq_Results/EUR_rs1061810_All_Genes.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

    # Obtain list of expressed repeats
    my.repeats <- Repeat_family.gs[Repeat_family.gs$gene %in% rownames(DESeq.res), 'gene']

    # Add the log2FC to the df made earlier
    DESeq.Condition3 <- DESeq.Condition3[my.repeats, ]
    DESeq.Condition3$Change <- DESeq.res[my.repeats, 'log2FoldChange']

# Load SNV4 (Condition 4) DESeq results
DESeq.res <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/1_DESeq_Genotypes/DESeq_Results/EUR_rs112581165_All_Genes.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

    # Obtain list of expressed repeats
    my.repeats <- Repeat_family.gs[Repeat_family.gs$gene %in% rownames(DESeq.res), 'gene']

    # Add the log2FC to the df made earlier
    DESeq.Condition4 <- DESeq.Condition4[my.repeats, ]
    DESeq.Condition4$Change <- DESeq.res[my.repeats, 'log2FoldChange']

# Load SNV5 (Condition 5) DESeq results
DESeq.res <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/1_DESeq_Genotypes/DESeq_Results/EUR_rs72691418_All_Genes.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

    # Obtain list of expressed repeats
    my.repeats <- Repeat_family.gs[Repeat_family.gs$gene %in% rownames(DESeq.res), 'gene']

    # Add the log2FC to the df made earlier
    DESeq.Condition5 <- DESeq.Condition5[my.repeats, ]
    DESeq.Condition5$Change <- DESeq.res[my.repeats, 'log2FoldChange']

# Organize family names alphabetically
DESeq.Condition1 <- DESeq.Condition1[order(DESeq.Condition1$gs_name), ]
DESeq.Condition2 <- DESeq.Condition2[order(DESeq.Condition2$gs_name), ]
DESeq.Condition3 <- DESeq.Condition3[order(DESeq.Condition3$gs_name), ]
DESeq.Condition4 <- DESeq.Condition4[order(DESeq.Condition4$gs_name), ]
DESeq.Condition5 <- DESeq.Condition5[order(DESeq.Condition5$gs_name), ]

# Lock in factors
DESeq.Condition1$gs_name <- factor(DESeq.Condition1$gs_name , levels = GSEA.families)
DESeq.Condition2$gs_name <- factor(DESeq.Condition2$gs_name , levels = GSEA.families)
DESeq.Condition3$gs_name <- factor(DESeq.Condition3$gs_name , levels = GSEA.families)
DESeq.Condition4$gs_name <- factor(DESeq.Condition4$gs_name , levels = GSEA.families)
DESeq.Condition5$gs_name <- factor(DESeq.Condition5$gs_name , levels = GSEA.families)
    
    
    


    

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
stats.condition3 <- calc_pval(input.DESeq.Condition = DESeq.Condition3)
stats.condition4 <- calc_pval(input.DESeq.Condition = DESeq.Condition4)
stats.condition5 <- calc_pval(input.DESeq.Condition = DESeq.Condition5)







# GENERATE PLOTS
my.label.yval <- -0.35

# IL16/STARD5 SNV
pdf(paste(my.output.dir, "Boxplot_Repeat_log2FC_EUR_SNV_IL16", ".pdf", sep=""), width = 3.25, height = 4)

    # Remove white space from header (3rd value in vector)
    par(mar=c(5,4,1,2)+0.1)

    # make boxplots for genotype segregated gene expression
    boxplot(Change ~ gs_name, 
            data = DESeq.Condition1,
            ylim = c(-0.4, 0.4),
            outline = F, 
            #col = boxplot.col,
            main = c('rs11635336 Effects (IL16/STARD5)'),
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
    #text(1.5, 9000, signif(OF.YF$p.value, 3), cex = 0.75)
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
    


# HLA SNV
my.label.yval <- -0.25

pdf(paste(my.output.dir, "Boxplot_Repeat_log2FC_EUR_SNV_HLA", ".pdf", sep=""), width = 3.25, height = 4)

    # Remove white space from header (3rd value in vector)
    par(mar=c(5,4,1,2)+0.1)

    # make boxplots for genotype segregated gene expression
    boxplot(Change ~ gs_name, 
            data = DESeq.Condition2,
            ylim = c(-0.3, 0.3),
            outline = F, 
            #col = boxplot.col,
            main = c('rs9271894 Effects (HLA)'),
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
    #text(1.5, 9000, signif(OF.YF$p.value, 3), cex = 0.75)
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
    


# HSD17B12 SNV
my.label.yval <- -0.25

pdf(paste(my.output.dir, "Boxplot_Repeat_log2FC_EUR_SNV_HSD17B12", ".pdf", sep=""), width = 3.25, height = 4)

    # Remove white space from header (3rd value in vector)
    par(mar=c(5,4,1,2)+0.1)

    # make boxplots for genotype segregated gene expression
    boxplot(Change ~ gs_name, 
            data = DESeq.Condition3,
            ylim = c(-0.3, 0.3),
            outline = F, 
            #col = boxplot.col,
            main = c('rs1061810 Effects (HSD17B12)'),
            cex.main = 0.70,
            ylab = c('log2(Fold Change)'),
            xlab = c(''),
            cex.lab = 0.7,
            cex.axis = 0.7,
            xaxt = 'n')
    
    # overlay the points
    stripchart(Change ~ gs_name, data = DESeq.Condition3, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 0.3, col = 'red') 
    
    # Add line at y = 0
    abline(h = 0, lty = 5, col = 'black')
    
    # Specify axis labels
    axis(1, at = 1:10, labels = GSEA.families, las = 2, cex.axis = 0.70)

    # Add pvalues
    #text(1.5, 9000, signif(OF.YF$p.value, 3), cex = 0.75)
    text(1, my.label.yval, stats.condition3[1, 'FDR_short'], cex = 0.7, srt = 90)
    text(2, my.label.yval, stats.condition3[2, 'FDR_short'], cex = 0.7, srt = 90)
    text(3, my.label.yval, stats.condition3[3, 'FDR_short'], cex = 0.7, srt = 90)
    text(4, my.label.yval, stats.condition3[4, 'FDR_short'], cex = 0.7, srt = 90)
    text(5, my.label.yval, stats.condition3[5, 'FDR_short'], cex = 0.7, srt = 90)
    text(6, my.label.yval, stats.condition3[6, 'FDR_short'], cex = 0.7, srt = 90)
    text(7, my.label.yval, stats.condition3[7, 'FDR_short'], cex = 0.7, srt = 90)
    text(8, my.label.yval, stats.condition3[8, 'FDR_short'], cex = 0.7, srt = 90)
    text(9, my.label.yval, stats.condition3[9, 'FDR_short'], cex = 0.7, srt = 90)
    text(10, my.label.yval, stats.condition3[10, 'FDR_short'], cex = 0.7, srt = 90)

    
# End pdf
dev.off()        
    


# SNV4
my.label.yval <- -0.55

pdf(paste(my.output.dir, "Boxplot_Repeat_log2FC_EUR_SNV4", ".pdf", sep=""), width = 3.25, height = 4)

    # Remove white space from header (3rd value in vector)
    par(mar=c(5,4,1,2)+0.1)

    # make boxplots for genotype segregated gene expression
    boxplot(Change ~ gs_name, 
            data = DESeq.Condition4,
            ylim = c(-0.6, 0.2),
            outline = F, 
            #col = boxplot.col,
            main = c('rs112581165 Effects (SNV4)'),
            cex.main = 0.70,
            ylab = c('log2(Fold Change)'),
            xlab = c(''),
            cex.lab = 0.7,
            cex.axis = 0.7,
            xaxt = 'n')
    
    # overlay the points
    stripchart(Change ~ gs_name, data = DESeq.Condition4, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 0.3, col = 'red') 
    
    # Add line at y = 0
    abline(h = 0, lty = 5, col = 'black')
    
    # Specify axis labels
    axis(1, at = 1:10, labels = GSEA.families, las = 2, cex.axis = 0.70)

    # Add pvalues
    #text(1.5, 9000, signif(OF.YF$p.value, 3), cex = 0.75)
    text(1, my.label.yval, stats.condition4[1, 'FDR_short'], cex = 0.7, srt = 90)
    text(2, my.label.yval, stats.condition4[2, 'FDR_short'], cex = 0.7, srt = 90)
    text(3, my.label.yval, stats.condition4[3, 'FDR_short'], cex = 0.7, srt = 90)
    text(4, my.label.yval, stats.condition4[4, 'FDR_short'], cex = 0.7, srt = 90)
    text(5, my.label.yval, stats.condition4[5, 'FDR_short'], cex = 0.7, srt = 90)
    text(6, my.label.yval, stats.condition4[6, 'FDR_short'], cex = 0.7, srt = 90)
    text(7, my.label.yval, stats.condition4[7, 'FDR_short'], cex = 0.7, srt = 90)
    text(8, my.label.yval, stats.condition4[8, 'FDR_short'], cex = 0.7, srt = 90)
    text(9, my.label.yval, stats.condition4[9, 'FDR_short'], cex = 0.7, srt = 90)
    text(10, my.label.yval, stats.condition4[10, 'FDR_short'], cex = 0.7, srt = 90)

    
# End pdf
dev.off()            
    


# SNV5
my.label.yval <- -0.50

pdf(paste(my.output.dir, "Boxplot_Repeat_log2FC_EUR_SNV5", ".pdf", sep=""), width = 3.25, height = 4)

    # Remove white space from header (3rd value in vector)
    par(mar=c(5,4,1,2)+0.1)

    # make boxplots for genotype segregated gene expression
    boxplot(Change ~ gs_name, 
            data = DESeq.Condition5,
            ylim = c(-0.6, 1.2),
            outline = F, 
            #col = boxplot.col,
            main = c('rs72691418 Effects (SNV5)'),
            cex.main = 0.70,
            ylab = c('log2(Fold Change)'),
            xlab = c(''),
            cex.lab = 0.7,
            cex.axis = 0.7,
            xaxt = 'n')
    
    # overlay the points
    stripchart(Change ~ gs_name, data = DESeq.Condition5, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 0.3, col = 'red') 
    
    # Add line at y = 0
    abline(h = 0, lty = 5, col = 'black')
    
    # Specify axis labels
    axis(1, at = 1:10, labels = GSEA.families, las = 2, cex.axis = 0.70)

    # Add pvalues
    #text(1.5, 9000, signif(OF.YF$p.value, 3), cex = 0.75)
    text(1, my.label.yval, stats.condition5[1, 'FDR_short'], cex = 0.7, srt = 90)
    text(2, my.label.yval, stats.condition5[2, 'FDR_short'], cex = 0.7, srt = 90)
    text(3, my.label.yval, stats.condition5[3, 'FDR_short'], cex = 0.7, srt = 90)
    text(4, my.label.yval, stats.condition5[4, 'FDR_short'], cex = 0.7, srt = 90)
    text(5, my.label.yval, stats.condition5[5, 'FDR_short'], cex = 0.7, srt = 90)
    text(6, my.label.yval, stats.condition5[6, 'FDR_short'], cex = 0.7, srt = 90)
    text(7, my.label.yval, stats.condition5[7, 'FDR_short'], cex = 0.7, srt = 90)
    text(8, my.label.yval, stats.condition5[8, 'FDR_short'], cex = 0.7, srt = 90)
    text(9, my.label.yval, stats.condition5[9, 'FDR_short'], cex = 0.7, srt = 90)
    text(10, my.label.yval, stats.condition5[10, 'FDR_short'], cex = 0.7, srt = 90)

    
# End pdf
dev.off()             
    
    
    
    




    
# Clean the environment
rm(list=ls())
    
    
