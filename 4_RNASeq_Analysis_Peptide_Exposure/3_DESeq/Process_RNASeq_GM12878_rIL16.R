# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/3_DESeq/Process_RNASeq_functions.R")

# Load libraries
library(DESeq2) # For differential expression
library(BiocParallel) # Used by DESeq2 for parallelization
  register(MulticoreParam(6)) # Use six cores
library(sva) # Correcting for batch effects
library(limma) # Correcting for batch effects
library(pheatmap) # For gene expression heatmaps
library(biomaRt) # To map Ensembl IDs to gene symbols
  

    
    
    
# Section 1: FILTER LOW EXPRESSION GENES ------------------------------------------------------------------------------------------------------------------------
 


  
# Define the output directory
dir.output <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/3_DESeq/Processed_counts_GM12878_rIL16/'

# Define the directory with the counts tables and collect individual file paths in a list.
count_Table.dir <- c("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/2_Read_counting/Counts_GM12878_rIL16/")
cntTable_list <- list.files(count_Table.dir, "\\.cntTable$", recursive=FALSE, full.names=TRUE)

# Run function to combine the individual counts tables into one df.
gene_TE_counts_raw <- aggregate_counts(cntTable_complete.list = cntTable_list)

# Remove useless info. from sample names
colnames(gene_TE_counts_raw) <- sub("Aligned.sortedByCoord.out.bam*", "", colnames(gene_TE_counts_raw))
colnames(gene_TE_counts_raw) <- sub("*fastp_", "", colnames(gene_TE_counts_raw))

# Reorder columns so controls are first
gene_TE_counts_raw <- gene_TE_counts_raw[, c('gene.TE', 
                                               'JB17', 'JB20', 'JB23',
                                               'JB18', 'JB21', 'JB24',
                                               'JB19', 'JB22', 'JB25')]
      
## Assign more specific sample names
colnames(gene_TE_counts_raw) <- c('gene.TE', 
                                  'T0_A', 'T0_B', 'T0_C',
                                  'T24_A', 'T24_B', 'T24_C', 
                                  'T48_A', 'T48_B', 'T48_C')

# Remove non-gene/TE features from expression table
misc.feature.indices <- which(!grepl('EBV', gene_TE_counts_raw$gene.TE) & !grepl('pcDNA', gene_TE_counts_raw$gene.TE))
gene_TE_counts_raw <- gene_TE_counts_raw[misc.feature.indices, ]

# Run function to remove gene/transcript version info, combine transcripts from the same gene, and filter lowly expressed genes *IMPORTANT: PSEUDOAUTOSOMAL GENES (PAR_Y) ARE COMBINED
counts_filtered <- cleanup_and_filter_counts(count_data = gene_TE_counts_raw, filter = TRUE, min_counts_per_sample = 10, fraction_of_samples = 3/9)

    # Save counts files.
    write.table(counts_filtered, file = paste(dir.output, "All_counts_filtered", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)

    




# Section 2: DESeq + MDS/PCA ------------------------------------------------------------------------------------------------------------------------




# Define the output directories
counts.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/3_DESeq/Processed_counts_GM12878_rIL16/'
DESeq.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/3_DESeq/DESeq_Results_GM12878_rIL16/'
MDS.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/3_DESeq/MDS_PCA_GM12878_rIL16/'



    
    
# DIFFERENTIAL EXPRESSION ANALYSIS (DESEQ2)
    
  # Define the counts matrix
  gene_TE_counts <- counts_filtered

  # General parameters
  padj_limit <- 0.05 # Define alpha: the significance cutoff used for optimizing the independent filtering (by default 0.1).   If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
  
  # Define batch effects and covariates (!!!!!! MAKE SURE VALUES MATCH IF COLUMN ORDER IS CHANGED !!!!!!)
  treatment.group <- c(rep('T0', 3), rep('T24', 3), rep('T48', 3))

  # Collect covariates in a dataframe
  SampleInfo <- data.frame(row.names = colnames(gene_TE_counts),
                               Treatment = treatment.group)
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = gene_TE_counts,
                                    colData = SampleInfo,
                                    design = ~ Treatment) 
  
  
  # run DESeq2
  dds <- DESeq(dds, parallel = TRUE)

  # Get VST data (includes size factor normalization)
  normalized_counts <- assay(varianceStabilizingTransformation(dds, blind = FALSE))

       # Save VST data
       write.table(normalized_counts, file = paste(counts.output.dir, 'All_counts_filtered_VST', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)

  
  #Extract DESeq results. 
  T24_vs_T0 <- results(dds, contrast = c("Treatment", "T24", "T0"), alpha = padj_limit, independentFiltering = TRUE)
  T48_vs_T0 <- results(dds, contrast = c("Treatment", "T48", "T0"), alpha = padj_limit, independentFiltering = TRUE)
  T48_vs_T24 <- results(dds, contrast = c("Treatment", "T48", "T24"), alpha = padj_limit, independentFiltering = TRUE)

  # DESeq Stats
  summary(T24_vs_T0, alpha = padj_limit)
  summary(T48_vs_T0, alpha = padj_limit)
  summary(T48_vs_T24, alpha = padj_limit)

  # Extract significant results and save.
  T24_vs_T0.sig <- Extract_DESeq_stats(DESeq_results = T24_vs_T0, padj_limit = padj_limit, organism = 'hs', 
                                                output.dir = DESeq.output.dir,
                                                output_file_prefix_all = 'rIL16_T24_vs_T0_All_Genes',
                                                output_file_prefix_sig = 'rIL16_T24_vs_T0_FDR5')
  
  T48_vs_T0.sig <- Extract_DESeq_stats(DESeq_results = T48_vs_T0, padj_limit = padj_limit, organism = 'hs', 
                                                output.dir = DESeq.output.dir,
                                                output_file_prefix_all = 'rIL16_T48_vs_T0_All_Genes',
                                                output_file_prefix_sig = 'rIL16_T48_vs_T0_FDR5')
  
  T48_vs_T24.sig <- Extract_DESeq_stats(DESeq_results = T48_vs_T24, padj_limit = padj_limit, organism = 'hs', 
                                                output.dir = DESeq.output.dir,
                                                output_file_prefix_all = 'rIL16_T48_vs_T24_All_Genes',
                                                output_file_prefix_sig = 'rIL16_T48_vs_T24_FDR5')

    
    

      
# CHECK FOR BATCHES BY MDS 

  # Define the expression data to analyze
  gene_TE_counts_normalized <- normalized_counts

  # Generate Gene-Only and TE-Only counts to assess their ability to stratify samples
  Gene_counts_normalized <- gene_TE_counts_normalized[grepl('ENSG', rownames(gene_TE_counts_normalized)), ]
  TE_counts_normalized <- gene_TE_counts_normalized[!grepl('ENSG', rownames(gene_TE_counts_normalized)), ]
  L1_expr <- gene_TE_counts_normalized[grepl(':L1:', rownames(gene_TE_counts_normalized)), ]

  # Run MDS
  mds.genes <- cmdscale(1-cor(Gene_counts_normalized, method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  mds.TE <- cmdscale(1-cor(TE_counts_normalized, method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  mds.L1 <- cmdscale(1-cor(L1_expr, method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)

  # Define MDS point parameters
  mds.colors <- c(rep('black', 3), rep('red1', 3), rep('red3', 3))
  mds.point.cex <- 1
  mds.label.cex <- 1
  mds.axis.cex <- 1
  mds.main.cex <- 1

  # Generate and save plot
  pdf(paste(MDS.output.dir, "Plot_MDS_rIL16", ".pdf", sep=""), width = 10, height = 10)
  par(mfrow=c(3,3))


        # Gene data
        plot(x = mds.genes[, 1],
             y = mds.genes[, 2],
             pch = 16,
             col = mds.colors,
             xlab = "MDS dimension 1",
             ylab = "MDS dimension 2",
             main = "Recombinant IL-16 Exposure, Gene Expression MDS",
             #xlim = c(-0.015, 0.015),
             #ylim = c(-0.010, 0.010),
             cex.main = mds.main.cex,
             cex.axis = mds.axis.cex,
             cex.lab = mds.label.cex,
             cex = mds.point.cex)

        legend("topleft", c("T0", "T24", "T48"), col = c('black', 'red1', 'red3'), pch = 16, pt.cex = 1, cex = 1, lty = 0, bty = 'n')
        
        # TE data
        plot(x = mds.TE[, 1],
             y = mds.TE[, 2],
             pch = 16,
             col = mds.colors,
             xlab = "MDS dimension 1",
             ylab = "MDS dimension 2",
             main = "Recombinant IL-16 Exposure, Repeat Expression MDS",
             #xlim = c(-0.02, 0.02),
             #ylim = c(-0.02, 0.02),
             cex.main = mds.main.cex,
             cex.axis = mds.axis.cex,
             cex.lab = mds.label.cex,
             cex = mds.point.cex)


        # L1 data
        plot(x = mds.L1[, 1],
             y = mds.L1[, 2],
             pch = 16,
             col = mds.colors,
             xlab = "MDS dimension 1",
             ylab = "MDS dimension 2",
             main = "Recombinant IL-16 Exposure, L1 Expression MDS",
             #xlim = c(-0.02, 0.02),
             #ylim = c(-0.02, 0.02),
             cex.main = mds.main.cex,
             cex.axis = mds.axis.cex,
             cex.lab = mds.label.cex,
             cex = mds.point.cex)


  dev.off()

    
  
  
  

  
  
  
  
  
  
  
# CHECK GROUPS BY PCA

  # Check that var for each gene is > 0
  var.genes <- apply(Gene_counts_normalized, 1, var) > 0
  var.TEs <- apply(TE_counts_normalized, 1, var) > 0
  var.L1 <- apply(L1_expr, 1, var) > 0

  # Perform PCA using genes with variation
  PCA.genes <- prcomp(t(Gene_counts_normalized[var.genes, ]), scale = TRUE)
  PCA.TEs <- prcomp(t(TE_counts_normalized[var.TEs, ]), scale = TRUE)
  PCA.L1 <- prcomp(t(L1_expr[var.L1, ]), scale = TRUE)

  # Extract pca stats
  summary.PCA.genes <- summary(PCA.genes)
  summary.PCA.TEs <- summary(PCA.TEs)
  summary.PCA.L1 <- summary(PCA.L1)

  # Define PCA point parameters
  pca.colors <- c(rep('black', 3), rep('red1', 3), rep('red3', 3))
  

  # Generate and save plot
  pdf(paste(MDS.output.dir, "Plot_PCA_rIL16", ".pdf", sep=""), height = 10, width = 10)
  par(mfrow=c(3,3)) # Added to produce images of similar shape across analyses


      # PCA for Genes ONLY
      plot(PCA.genes$x[, 1],
           PCA.genes$x[, 2],
           col = pca.colors,
           main = "Recombinant IL-16 Exposure, Gene Expression PCA",
           xlab = paste('PC1 (', round(100*summary.PCA.genes$importance[,1][2],1),"%)", sep=""),
           ylab = paste('PC2 (', round(100*summary.PCA.genes$importance[,2][2],1),"%)", sep=""),
           #xlim = c(-200, 200),
           #ylim = c(-175, 175),
           cex.main = mds.main.cex,
           cex.lab = mds.label.cex,
           cex.axis = mds.axis.cex,
           pch = 16,
           cex = mds.point.cex)

      legend("bottomright", c("T0", "T24", "T48"), col = c('black', 'red1', 'red3'), pch = 16, pt.cex = 1, cex = 1, lty = 0, bty = 'n')

      # PCA for TEs ONLY
      plot(PCA.TEs$x[, 1],
           PCA.TEs$x[, 2],
           col = pca.colors,
           main = "Recombinant IL-16 Exposure, Repeat Expression PCA",
           xlab = paste('PC1 (', round(100*summary.PCA.TEs$importance[,1][2],1),"%)", sep=""),
           ylab = paste('PC2 (', round(100*summary.PCA.TEs$importance[,2][2],1),"%)", sep=""),
           #xlim = c(-100, 100),
           #ylim = c(-50, 50),
           cex.main = mds.main.cex,
           cex.lab = mds.label.cex,
           cex.axis = mds.axis.cex,
           pch = 16,
           cex = mds.point.cex)
      
      # PCA for L1s ONLY
      plot(PCA.L1$x[, 1],
           PCA.L1$x[, 2],
           col = pca.colors,
           main = "Recombinant IL-16 Exposure, L1 Expression PCA",
           xlab = paste('PC1 (', round(100*summary.PCA.TEs$importance[,1][2],1),"%)", sep=""),
           ylab = paste('PC2 (', round(100*summary.PCA.TEs$importance[,2][2],1),"%)", sep=""),
           #xlim = c(-100, 100),
           #ylim = c(-50, 50),
           cex.main = mds.main.cex,
           cex.lab = mds.label.cex,
           cex.axis = mds.axis.cex,
           pch = 16,
           cex = mds.point.cex)

  dev.off()
  
  
  
  


# Section 3: COMPARISON OF GENE EXPRESSION SCATTERPLOTS (GM12878 IL16 OE vs rIL16 T24/T48)  -------------------------------------------------



# directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/3_DESeq/DESeq_Results_GM12878_rIL16/'




# PREPARE DATA FOR T24 COMPARISON

# Load unfiltered data
ALL.IL16.OE <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/3_DESeq/DESeq_Results_GM12878_OE/All_Genes_IL16.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
ALL.IL16.T24 <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/3_DESeq/DESeq_Results_GM12878_rIL16/rIL16_T24_vs_T0_All_Genes.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Find shared gene sets
shared.genes <- which(rownames(ALL.IL16.OE) %in% rownames(ALL.IL16.T24))
shared.genes <- rownames(ALL.IL16.OE)[shared.genes]

# Keep data for common gene sets
ALL.IL16.OE  <- ALL.IL16.OE[shared.genes, ]
ALL.IL16.T24 <- ALL.IL16.T24[shared.genes, ]
          

      
                      
# GENERATE SCATTERPLOTS FOR T24 COMPARISON
    
# Spearman correlation between common gene sets 
gene.cor <- cor.test(x = rank(ALL.IL16.OE$log2FoldChange), y = rank(ALL.IL16.T24$log2FoldChange), method = 'spearman', exact = FALSE)

# Define cor and pval values to include in the legend
gene.rho  <- format(round(gene.cor[["estimate"]][["rho"]], 2), nsmall = 2)
gene.pval  <- formatC(gene.cor[["p.value"]], format = 'E', digits = 2)

# Generate the GO scatterplot
pdf(paste(my.output.dir, "Scatterplot_Individual_Genes_IL16_OE_vs_Peptide_T24", ".pdf", sep=""), height = 4, width = 4)

      # Smoothscatter plot 
      smoothScatter(ALL.IL16.OE$log2FoldChange, 
                    ALL.IL16.T24$log2FoldChange,
                    colramp = colorRampPalette(c("white", 'red2', 'red4')),
                    nrpoints = 0,
                    xlim = c(-2.5, 2.5),
                    ylim = c(-2.5, 2.5),
                    xlab = c('log2FC IL16 OE'),
                    ylab = c('log2FC rhIL16 T24'),
                    main = c('DEG Correlation'))
      
      # Add line for x=y
      abline(a = 0, b = 1, col = 'darkgrey', lty = 2)
      
      # add line for y = 0
      abline(h = 0, lty = 2)
      
      # add line for x = 0
      abline(v = 0, lty = 2)
      
      #add legend
      legend("topleft", 
             c(paste('rho = ', gene.rho), paste('p = ', gene.pval)),
             cex = 0.75)

dev.off()




# PREPARE DATA FOR T48 COMPARISON

# Load unfiltered data
ALL.IL16.OE <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/3_DESeq/DESeq_Results_GM12878_OE/All_Genes_IL16.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
ALL.IL16.T48 <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/3_DESeq/DESeq_Results_GM12878_rIL16/rIL16_T48_vs_T0_All_Genes.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Find shared gene sets
shared.genes <- which(rownames(ALL.IL16.OE) %in% rownames(ALL.IL16.T48))
shared.genes <- rownames(ALL.IL16.OE)[shared.genes]

# Keep data for common gene sets
ALL.IL16.OE  <- ALL.IL16.OE[shared.genes, ]
ALL.IL16.T48 <- ALL.IL16.T48[shared.genes, ]
          
       
      
                      
# GENERATE SCATTERPLOTS FOR T48 COMPARISON
    
# Spearman correlation between common gene sets 
gene.cor <- cor.test(x = rank(ALL.IL16.OE$log2FoldChange), y = rank(ALL.IL16.T48$log2FoldChange), method = 'spearman', exact = FALSE)

# Define cor and pval values to include in the legend
gene.rho  <- format(round(gene.cor[["estimate"]][["rho"]], 2), nsmall = 2)
gene.pval  <- formatC(gene.cor[["p.value"]], format = 'E', digits = 2)

# Generate the GO scatterplot
pdf(paste(my.output.dir, "Scatterplot_Individual_Genes_IL16_OE_vs_Peptide_T48", ".pdf", sep=""), height = 4, width = 4)

      # Smoothscatter plot 
      smoothScatter(ALL.IL16.OE$log2FoldChange, 
                    ALL.IL16.T48$log2FoldChange,
                    colramp = colorRampPalette(c("white", 'red2', 'red4')),
                    nrpoints = 0,
                    xlim = c(-2.5, 2.5),
                    ylim = c(-2.5, 2.5),
                    xlab = c('log2FC IL16 OE'),
                    ylab = c('log2FC rhIL16 T48'),
                    main = c('DEG Correlation'))
      
      # Add line for x=y
      abline(a = 0, b = 1, col = 'darkgrey', lty = 2)
      
      # add line for y = 0
      abline(h = 0, lty = 2)
      
      # add line for x = 0
      abline(v = 0, lty = 2)
      
      #add legend
      legend("topleft", 
             c(paste('rho = ', gene.rho), paste('p = ', gene.pval)),
             cex = 0.75)

dev.off()


  



# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/3_DESeq/Session_Info/'
    
    
sink(file = paste(dir.session_info,"Session_Info_GM12878_rIL16.txt", sep =""))
sessionInfo()
sink()      
    



    
# Clean the environment
rm(list=ls())
