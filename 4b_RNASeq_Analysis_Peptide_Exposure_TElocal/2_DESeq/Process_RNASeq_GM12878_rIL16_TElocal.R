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
dir.output <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/2_DESeq/Processed_counts_GM12878_rIL16/'

# Define the directory with the counts tables and collect individual file paths in a list.
count_Table.dir <- c("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/1_Read_counting/TElocal_Counts_GM12878_rIL16/")
cntTable_list <- list.files(count_Table.dir, "\\.cntTable$", recursive=FALSE, full.names=TRUE)

# Run function to combine the individual counts tables into one df.
gene_TE_counts_raw <- aggregate_counts(cntTable_complete.list = cntTable_list)

# Remove genes/TEs with all zeros
gene_TE_counts_raw <- gene_TE_counts_raw[which(rowSums(gene_TE_counts_raw[, -c(1)]) > 0), ]

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

# Load repeat subsets
repeats.intergenic.nearby <- read.csv("/Users/juanb/Documents/Bioinformatic_Tools/STAR_genome_indices/TElocal_rmsk_files/GRCh38_GENCODE_rmsk_TE_subsets_vs_GENCODE33/Repeats_intergenic_within_5kb.gtf", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)
repeats.intergenic.distal <- read.csv("/Users/juanb/Documents/Bioinformatic_Tools/STAR_genome_indices/TElocal_rmsk_files/GRCh38_GENCODE_rmsk_TE_subsets_vs_GENCODE33/Repeats_intergenic_farther_than_5kb.gtf", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)
repeats.intronic <- read.csv("/Users/juanb/Documents/Bioinformatic_Tools/STAR_genome_indices/TElocal_rmsk_files/GRCh38_GENCODE_rmsk_TE_subsets_vs_GENCODE33/Repeats_intronic.gtf", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)
repeats.exonic <- read.csv("/Users/juanb/Documents/Bioinformatic_Tools/STAR_genome_indices/TElocal_rmsk_files/GRCh38_GENCODE_rmsk_TE_subsets_vs_GENCODE33/Repeats_exonic.gtf", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)

    # Assign repeat_name to rownames
    rownames(repeats.intergenic.nearby) <- repeats.intergenic.nearby$V9
    rownames(repeats.intergenic.distal) <- repeats.intergenic.distal$V9
    rownames(repeats.intronic) <- repeats.intronic$V9
    rownames(repeats.exonic) <- repeats.exonic$V9
  
# Combine TE loci by various parameters (intergenic/exonic/intronic)
gene_TE_counts_raw_split <- prefilter_TE_aggregation(gene_TE_counts = gene_TE_counts_raw, 
                                                     aggregate_all = FALSE, 
                                                     intergenic_nearby_repeats = repeats.intergenic.nearby, 
                                                     intergenic_distal_repeats = repeats.intergenic.distal, 
                                                     intronic_repeats = repeats.intronic, 
                                                     exonic_repeats = repeats.exonic)

# Run function to remove gene/transcript version info, combine transcripts from the same gene, and filter lowly expressed genes *IMPORTANT: PSEUDOAUTOSOMAL GENES (PAR_Y) ARE COMBINED
counts_filtered <- cleanup_and_filter_counts(count_data = gene_TE_counts_raw_split, filter = TRUE, min_counts_per_sample = 10, fraction_of_samples = 3/9)

    # Save counts files.
    write.table(counts_filtered, file = paste(dir.output, "TElocal_counts_filtered", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)




# Section 2: DESeq + MDS/PCA ------------------------------------------------------------------------------------------------------------------------




# Define the output directories
counts.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/2_DESeq/Processed_counts_GM12878_rIL16/'
DESeq.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/2_DESeq/DESeq_Results_GM12878_rIL16/'
MDS.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/2_DESeq/MDS_PCA_GM12878_rIL16/'    



    
    
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
       write.table(normalized_counts, file = paste(counts.output.dir, 'TElocal_counts_filtered_VST', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)

  
  # Extract DESeq results. 
  T24_vs_T0 <- results(dds, contrast = c("Treatment", "T24", "T0"), alpha = padj_limit, independentFiltering = TRUE)
  T48_vs_T0 <- results(dds, contrast = c("Treatment", "T48", "T0"), alpha = padj_limit, independentFiltering = TRUE)

  # DESeq Stats
  summary(T24_vs_T0, alpha = padj_limit)
  summary(T48_vs_T0, alpha = padj_limit)

  # Extract significant results and save.
  T24_vs_T0.sig <- Extract_DESeq_stats(DESeq_results = T24_vs_T0, padj_limit = padj_limit, organism = 'hs', 
                                                output.dir = DESeq.output.dir,
                                                output_file_prefix_all = 'rIL16_T24_vs_T0_All_Genes',
                                                output_file_prefix_sig = 'rIL16_T24_vs_T0_FDR5')
  
  T48_vs_T0.sig <- Extract_DESeq_stats(DESeq_results = T48_vs_T0, padj_limit = padj_limit, organism = 'hs', 
                                                output.dir = DESeq.output.dir,
                                                output_file_prefix_all = 'rIL16_T48_vs_T0_All_Genes',
                                                output_file_prefix_sig = 'rIL16_T48_vs_T0_FDR5')


    
    

      
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
  
  
  
  


# Session Info  -------------------------------------------------


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/2_DESeq/Session_Info/'
    
    
sink(file = paste(dir.session_info,"Session_Info_GM12878_rIL16.txt", sep =""))
sessionInfo()
sink()      
    



    
# Clean the environment
rm(list=ls())
