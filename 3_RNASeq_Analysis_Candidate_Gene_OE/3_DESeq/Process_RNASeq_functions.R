aggregate_counts <- function(cntTable_complete.list) {
  
  # cntTable_complete.list should be a list of counts files
  
  
  
  # Read 1 count table and keep only the gene column. In the FOR loop below, extend by adding counts from each sample. 
  count_data <- read.csv(cntTable_complete.list[1], header=TRUE, sep="")
  count_data <- count_data[1]
  
  
  for (count_table_file in cntTable_complete.list) {
    
    # Load individual table temporarily. 
    count_data_temp <- read.csv(count_table_file, header=TRUE, sep="")
    
    # Check whether the GENE/TE names column matches that of the main count_data table.  
    column_1_match <- all.equal(count_data_temp[1], count_data[1])
    
    if (column_1_match == TRUE) {
      
      # Add the temporary counts data to the final counts object.
      count_data <- data.frame(count_data, count_data_temp[2])
      
      # Remove count_data_temp since it is not needed
      rm(count_data_temp)
      
    } else {
      
      count_data <- c('Count table gene names do not match in some samples')
      break
      
    } # End if else.
    
  } #End FOR loop
  
  
  
  # Output the final count data table.
  return(count_data)
  
  
  
  
} # END FUNCTION

cleanup_and_filter_counts <- function(count_data, filter, min_counts_per_sample, fraction_of_samples) {
  
  # FUNCTION INFO:
  # CLEANUP: REMOVE GENE/TRANSCRIPT VERSION AND COMBINE COUNTS FROM SIMILAR GENES/TRANSCRIPTS
  # FILTER: LOWLY EXPRESSED GENES
  
  

  # Assign name to columns with gene/TE names. Name needs to specifically be 'gene.TE' since that is referenced.
  names(count_data)[1] <- 'gene.TE'    

  
  # Remove transcript or ENSEMBL version info. Keeping this can raise issues with gene name mapping.
  count_data[[1]] <- sub("\\..*", "", count_data[[1]], fixed=FALSE)
  
  
  # Collect same transcripts into one row. NOTE, pseudoautosomal PAR_Y genes will be combined with genes in homologous regions.
  gene_counts_summed <- aggregate(count_data[,-c(1)],
                                  by = list(count_data$gene.TE),
                                  FUN = 'sum')
  
  # Update the counts matrix by moving the gene/TE column to rownames and deleting the original column
  rownames(gene_counts_summed) <- gene_counts_summed[, 1]
  gene_counts_summed <- gene_counts_summed[, c(-1)]
  
  
  # Remove lowly expressed genes, if desired. Filter by CPM (about 10 reads across all samples) to control for differences in library depth.
  if (filter == TRUE) {
    
    # Calculate the number of samples in the counts table
    num_of_samples <- ncol(count_data) - 1
    
    # Define the minimum number of samples that have to meet the CPM threshold
    min.samples <- round(fraction_of_samples * num_of_samples)
    
    # Calculate library sizes for each sample 
    library_sizes <- colSums(gene_counts_summed)
    
    # Calculate the median library size
    MedianLibSize <- median(library_sizes)
    
    # calculate the CPM cutoff for the input libraries (changes with the median library size)
    CPM.Cutoff <- min_counts_per_sample / MedianLibSize*1e6
    
    # Convert counts to CPMs
    CPM_table <- edgeR::cpm(gene_counts_summed, lib.size = library_sizes)
    
    # Make a filter function, where "min.samples" have to exceed "CPM.Cutoff"
    cpm.filter.func <- genefilter::kOverA(min.samples, CPM.Cutoff)
    
    # Bind the filtering function to filterfun
    flist <- genefilter::filterfun(cpm.filter.func)
    
    # Identify genes that pass the filtering threshold
    kept_genes <- genefilter::genefilter(CPM_table, flist)
    
    # Get counts for the genes that pass filtering
    gene_counts_summed_filtered <- gene_counts_summed[kept_genes, ]
    
    # Return the cleaned and filtered counts
    return(gene_counts_summed_filtered)
    
  } else {
    
    # Return the cleaned counts
    return(gene_counts_summed)
    
  }
  
  

  
} # END FUNCTION

Ensembl_to_symbol <- function(organism, Ensembl_genes){
  
  # THIS FUNCTION MAPS AN ENSEMBL GENE NAME TO ITS SYMBOL
  
  
  
  # Define organism specific parameters
  if (organism == 'hs') {
    organism_dataset <- "hsapiens_gene_ensembl"
    Ensembl_version <- 99
    
  } else if (organism == 'mm') {
    organism_dataset <- "mmusculus_gene_ensembl"
    Ensembl_version <- 99 # DOUBLE CHECK
  }
  
  
  
  # Retrieve ensembl database
  ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset = organism_dataset, host = 'https://www.ensembl.org', version = Ensembl_version, verbose = TRUE)
  
  # Map Ensembl names to gene symbols
  Mapped_gene_info <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol", "entrezgene_id"), 
                                     filters = c("ensembl_gene_id"), 
                                     values = Ensembl_genes,
                                     mart = ensembl,
                                     verbose = FALSE,
                                     uniqueRows = TRUE)  
  
  # If genes have multiple mappings, keep only the first. 
  Mapped_gene_info <- Mapped_gene_info[!duplicated(Mapped_gene_info$ensembl_gene_id), ]
  
  # Assign Ensembl names to the rownames
  rownames(Mapped_gene_info) <- Mapped_gene_info$ensembl_gene_id
  
  # Remove ensembl gene id column since they're in rownames
  Mapped_gene_info <- Mapped_gene_info[, -c(1)] 
  
  # Return table with the mapping info
  return(Mapped_gene_info)
  
  
} # END FUNCTION

Run_SVA <- function(filtered_counts, control_label, control_reps, treatment_label, treatment_reps, output_location){
  

  
  
  
  
  
# Make sure counts are in matrix format
filtered_counts <- as.matrix(filtered_counts)


# Define treatments (!!!!!! MAKE SURE VALUES MATCH IF COLUMN ORDER IS CHANGED !!!!!!)
treatment.group <- c(rep(control_label, control_reps), 
                   rep(treatment_label, treatment_reps)
                   )


# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = colnames(filtered_counts),
                         Treatment = treatment.group
                         )

# Set the full models (adjustment + variables of interest).
mod_full = model.matrix(~ Treatment, data = SampleInfo)

# Set null models (adjustment variables)
mod_null = model.matrix(~ 1, data = SampleInfo)

# estimate the # of SVs
n.sv.be = num.sv(dat = filtered_counts, mod = mod_full, method ="be", B = 50, seed = 12345) 

# apply SVAseq algorithm
my_svaseq = svaseq(dat = filtered_counts, mod = mod_full, mod0 = mod_null, n.sv = n.sv.be, constant = 0.1, B = 50)

# Remove batch effects with Limma
my.corrected.data <- removeBatchEffect(x = log2(filtered_counts + 0.1),
                                    batch = NULL,
                                    covariates = cbind(my_svaseq$sv), # !!!!!!! Change with SVs
                                    design = mod_full)

# Delog and round data.
my.corrected.data <- round(2^my.corrected.data-0.1)

    # Save batch corrected counts (to visualize downstream results)
    write.table(my.corrected.data, file = paste(output_location, 'All_counts_filtered_SVA_', treatment_label, '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)  




  
  

# Return the corrected data
return(my.corrected.data)
  
  
  
} # END FUNCTION

Known_Batch_Removal <- function(filtered_counts, control_label, control_reps, treatment_label, treatment_reps, output_location, batch_df){
  

  
  
  
  
  
# Make sure counts are in matrix format
filtered_counts <- as.matrix(filtered_counts)


# Define treatments (!!!!!! MAKE SURE VALUES MATCH IF COLUMN ORDER IS CHANGED !!!!!!)
treatment.group <- c(rep(control_label, control_reps), 
                   rep(treatment_label, treatment_reps)
                   )


# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = colnames(filtered_counts),
                         Treatment = treatment.group
                         )

# Set the full models (adjustment + variables of interest).
mod_full = model.matrix(~ Treatment, data = SampleInfo)

# Remove batch effects with Limma
my.corrected.data <- removeBatchEffect(x = log2(filtered_counts + 0.1),
                                    batch = NULL,
                                    covariates = cbind(batch_df), 
                                    design = mod_full)

# Delog and round data.
my.corrected.data <- round(2^my.corrected.data-0.1)

    # Save batch corrected counts (to visualize downstream results)
    write.table(my.corrected.data, file = paste(output_location, 'All_counts_filtered_batchesremoved_', treatment_label, '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)  




  
  

# Return the corrected data
return(my.corrected.data)
  
  
  
} # END FUNCTION

Run_DESeq <- function(filtered_counts, padj_limit, control_label, control_reps, treatment_label, treatment_reps, VST_output, DESeq_output){
  
# padj_limit defines alpha: the significance cutoff used for optimizing the independent filtering (by default 0.1).   If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.

  
  
  


# Define the treatments
treatment.group <- c(rep(control_label, control_reps), 
                     rep(treatment_label, treatment_reps)
                     )

# Collect covariates in a dataframe
SampleInfo <- data.frame(row.names = colnames(filtered_counts),
                         Treatment = treatment.group)


# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                                  colData = SampleInfo,
                                  design = ~ Treatment) 

# run DESeq2
dds <- DESeq(dds, parallel = TRUE)

# Get VST data (includes size factor normalization)
my.VST <- assay(varianceStabilizingTransformation(dds, blind = FALSE))

     # Save VST data
     write.table(my.VST, file = paste(VST_output, 'All_counts_filtered_', treatment_label, '_VST', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
     



# Extract DESeq results. 
DESeq.Effects <- results(dds, contrast = c("Treatment", treatment_label, control_label), alpha = padj_limit, independentFiltering = TRUE)


# DESeq Stats
summary(DESeq.Effects, alpha = padj_limit)



# Generate MA plots
pdf(paste(DESeq_output,"Plot_MA_DESeq_Results_", treatment_label, ".pdf", sep=""), height = 10, width = 10)

    DESeq2::plotMA(object = DESeq.Effects, alpha = padj_limit, main = 'MA Plot', colSig = 'orange', ylim = c(-4,4))

dev.off()





# Run function to map ensembl labels, extract significant results, and save them.
my.sig <- Extract_DESeq_stats(DESeq_results = DESeq.Effects, padj_limit = padj_limit, organism = 'hs', 
                              output.dir = DESeq_output,
                              output_file_prefix_all = paste('All_Genes_', treatment_label, sep = ''),
                              output_file_prefix_sig = paste('FDR_', padj_limit, '_Genes_', treatment_label, sep = ''))

  
  

# Return the VST data and significant results as a list

return(list(my.VST, my.sig))
  
  
  
} # END FUNCTION

Extract_DESeq_stats <- function(DESeq_results, padj_limit, organism, output.dir, output_file_prefix_all, output_file_prefix_sig){
  
  # FUNCTION NOTE: Extract statistics from DESeq Results object, remove genes with NA for padj, save padj filtered and unfiltered list. The unfiltered list can be used as gene background.
  
  
  
  # Check if output directory exists. If not, make it. 
  if ( dir.exists(output.dir) == FALSE) {
    dir.create(output.dir, recursive = TRUE)
  }
  
  
  
  # Extract statistics table from DESeq results
  DESeq_stats <- as.data.frame(DESeq_results@listData)
  
  # Assign rownames (genes) to statistics table
  rownames(DESeq_stats) <- DESeq_results@rownames
  
  # Remove rows with NAs
  DESeq_stats <- na.omit(DESeq_stats)
  
  # Getting alternative gene names for remaining genes
  alternative_names <- Ensembl_to_symbol(organism = organism, Ensembl_genes = rownames(DESeq_stats))
  
  # Make columns to hold alternative gene names
  DESeq_stats$external_gene_name <- NA 
  DESeq_stats$hgnc_symbol <- NA
  DESeq_stats$entrezgene_id <- NA
  
  # Assign alternative names
  DESeq_stats[rownames(alternative_names),c('external_gene_name', 'hgnc_symbol', 'entrezgene_id')] <- alternative_names[,]

  # Filter out DEGs
  DESeq_stats.sig <- DESeq_stats[DESeq_stats$padj < padj_limit, ]
  
  # Save stats for the full and significant gene names
  write.table(DESeq_stats.sig, file = paste(output.dir, output_file_prefix_sig, '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)
  write.table(DESeq_stats, file = paste(output.dir, output_file_prefix_all, '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)

  
  
  # Output the significant gene stats
  return(DESeq_stats.sig)
  
} # END FUNCTION

Plot_Expression <- function(output.dir, target_gene_Ensembl_name, target_gene_symbol, VST.data, Sig.DESeq.Res, control_label, control_n, treatment_label, treatment_n, ymin, ymax){
  
  
  
  
  
  
    # Make a vector with the group each sample belongs to
    sample.groups <- c(rep(control_label, control_n), rep(treatment_label, treatment_n))
    
    # Make a dataframe with the expression info
    expression <- data.frame(row.names = colnames(VST.data),
                             my.group = sample.groups,
                             target.expression = VST.data[target_gene_Ensembl_name, ])
    
    # Extract the padj value
    target.padj <- Sig.DESeq.Res[target_gene_Ensembl_name, 'padj']
  
    # Start PDF to plot other responses
    pdf(paste(output.dir, "Boxplot_Gene_Expression_", target_gene_symbol, ".pdf", sep=""), width = 4.4, height = 4.5)
    #par(mfrow = c(2, 2)) # 2x2 plotting matrix
        
            # make boxplots 
            boxplot(target.expression ~ my.group, 
                    data = expression,
                    ylim = c(ymin, ymax),
                    outline = F, 
                    main = paste(target_gene_symbol, ' Expression', sep = ''),
                    ylab = paste('Normalized log2(counts)', sep = ''),
                    xlab = c(''),
                    xaxt = 'n')
            
            # overlay the points
            stripchart(target.expression ~ my.group, data = expression, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = 'red') 
            
            # Specify axis labels
            axis(1, at = c(1, 2), labels = c(control_label, treatment_label))
            
            # Add pvalues
            text(1.5, ymax, paste('adj.p = ', signif(target.padj, 3), sep = ''), cex = 1)

    
    # End pdf
    dev.off()
  
  
  
    # No output needed
    return()
  
  
} # END FUNCTION

Candidate_gene_heatmaps <- function(output.dir, EnsemblIDs_to_plot, my.gene.symbols, VST.data, control_label, control_n, control_color, treatment_label, treatment_n, treatment_color){
  
  # NOTE: THIS FUNCTION RELIES ON PHEATMAP
  
  
  
  

  
  # Extract expression values for genes to plot
  my.gene.expression <- VST.data[EnsemblIDs_to_plot, ]
  
  # Change expression rownames to symbols
  rownames(my.gene.expression) <- my.gene.symbols
  
  # Define column label groups 
  plot_column_groups <- data.frame(row.names = colnames(VST.data), 
                                   Treatment = c(rep(control_label, control_n), rep('OE', treatment_n)))
      
  # Define group colors
  group_color_map <- c(control_color, treatment_color)
  names(group_color_map) <- c(control_label, 'OE')
  plot_column_colors <- list(Treatment = group_color_map)
  
  # Define pheatmap scale
  breaksList <- seq(-2, 2, by = 0.02)
      
  # Generate heatmap
  pdf(paste(output.dir, "Heatmap_Candidate_Gene_Coexpression_", treatment_label, ".pdf", sep=""), height = 3, width = 3)
  
  pheatmap(my.gene.expression, 
           cluster_rows = FALSE, 
           show_rownames = TRUE,
           fontsize_row = 9,
           fontsize = 9,
           cluster_cols = FALSE, 
           show_colnames = FALSE,
           scale = 'row',
           border_color = NA,
           annotation_col = plot_column_groups,
           annotation_colors = plot_column_colors,
           color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
           breaks = breaksList
           )

  dev.off()     
          
  
  
  
} # END FUNCTION

Run_MDS <-function(VST_expression, control_color, control_reps, treatment_label, treatment_color, treatment_reps, MDS_dir){
  
  
  
  
  
  
    
# Generate Gene-Only, TE-Only, LINE-only, and L1-only expression to assess their ability to stratify samples
Gene_expr <- VST_expression[grepl('ENSG', rownames(VST_expression)), ]
TE_expr <- VST_expression[!grepl('ENSG', rownames(VST_expression)) & !grepl('EBV_gene', rownames(VST_expression)) & !grepl('pcDNA_NeoR_gene', rownames(VST_expression)), ]
LINE_expr <- VST_expression[grepl(':LINE', rownames(VST_expression)), ]
L1_expr <- VST_expression[grepl(':L1:', rownames(VST_expression)), ]

# Run MDS
mds.genes <- cmdscale(1-cor(Gene_expr, method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
mds.TE <- cmdscale(1-cor(TE_expr, method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
mds.LINE <- cmdscale(1-cor(LINE_expr, method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
mds.L1 <- cmdscale(1-cor(L1_expr, method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)

# Define MDS point colors
mds.colors <- c(rep(control_color, control_reps), 
                rep(treatment_color, treatment_reps)
                )

# Define MDS point parameters
mds.point.cex <- 1
mds.label.cex <- 1
mds.axis.cex <- 1
mds.main.cex <- 1

# Generate and save plot
pdf(paste(MDS_dir, "Plot_MDS_", treatment_label, ".pdf", sep=""), width = 10, height = 10)
par(mfrow=c(3,3))

      # Gene data
      plot(x = mds.genes[, 1],
           y = mds.genes[, 2],
           pch = 16,
           col = mds.colors,
           xlab = "MDS dimension 1",
           ylab = "MDS dimension 2",
           main = paste(treatment_label, " OE, Gene Expression MDS", sep = ''),
           #xlim = c(-0.015, 0.015),
           #ylim = c(-0.010, 0.010),
           cex.main = mds.main.cex,
           cex.axis = mds.axis.cex,
           cex.lab = mds.label.cex,
           cex = mds.point.cex)

      # TE data
      plot(x = mds.TE[, 1],
           y = mds.TE[, 2],
           pch = 16,
           col = mds.colors,
           xlab = "MDS dimension 1",
           ylab = "MDS dimension 2",
           main = paste(treatment_label, " OE, Repeat Expression MDS", sep = ''),
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
           main = paste(treatment_label, " OE, L1 Expression MDS", sep = ''),
           #xlim = c(-0.02, 0.02),
           #ylim = c(-0.02, 0.02),
           cex.main = mds.main.cex,
           cex.axis = mds.axis.cex,
           cex.lab = mds.label.cex,
           cex = mds.point.cex)

dev.off()
  
  

# Return nothing
return(NULL)
  
  
  
} # END FUNCTION

Run_PCA <-function(VST_expression, control_color, control_reps, treatment_label, treatment_color, treatment_reps, MDS_dir){
  
  
  
  
  
  
    
# Generate Gene-Only, TE-Only, LINE-only, L1-only expression to assess their ability to stratify samples
Gene_expr <- VST_expression[grepl('ENSG', rownames(VST_expression)), ]
TE_expr <- VST_expression[!grepl('ENSG', rownames(VST_expression)) & !grepl('EBV_gene', rownames(VST_expression)) & !grepl('pcDNA_NeoR_gene', rownames(VST_expression)), ]
LINE_expr <- VST_expression[grepl(':LINE', rownames(VST_expression)), ]
L1_expr <- VST_expression[grepl(':L1:', rownames(VST_expression)), ]

# Check that var for each gene is > 0
var.genes <- apply(Gene_expr, 1, var) > 0
var.TEs <- apply(TE_expr, 1, var) > 0
var.LINE <- apply(LINE_expr, 1, var) > 0
var.L1 <- apply(L1_expr, 1, var) > 0

# Perform PCA using genes with variation
PCA.genes <- prcomp(t(Gene_expr[var.genes, ]), scale = TRUE)
PCA.TEs <- prcomp(t(TE_expr[var.TEs, ]), scale = TRUE)
PCA.LINE <- prcomp(t(LINE_expr[var.LINE, ]), scale = TRUE)
PCA.L1 <- prcomp(t(L1_expr[var.L1, ]), scale = TRUE)

# Extract pca stats
summary.PCA.genes <- summary(PCA.genes)
summary.PCA.TEs <- summary(PCA.TEs)
summary.PCA.LINE <- summary(PCA.LINE)
summary.PCA.L1 <- summary(PCA.L1)

# Define MDS point colors
pca.colors <- c(rep(control_color, control_reps), 
                rep(treatment_color, treatment_reps)
                )
mds.point.cex <- 1
mds.label.cex <- 1
mds.axis.cex <- 1
mds.main.cex <- 1


# Generate and save plot
pdf(paste(MDS_dir, "Plot_PCA_", treatment_label, ".pdf", sep=""), width = 10, height = 10)
par(mfrow=c(3,3)) # Added to produce images of similar shape across analyses

    # PCA for Genes ONLY
    plot(PCA.genes$x[, 1],
         PCA.genes$x[, 2],
         col = pca.colors,
         main = paste(treatment_label, " OE, Gene Expression PCA", sep = ''),
         xlab = paste('PC1 (', round(100*summary.PCA.genes$importance[,1][2],1),"%)", sep=""),
         ylab = paste('PC2 (', round(100*summary.PCA.genes$importance[,2][2],1),"%)", sep=""),
         #xlim = c(-200, 200),
         #ylim = c(-90, 90),
         cex.main = mds.main.cex,
         cex.lab = mds.label.cex,
         cex.axis = mds.axis.cex,
         pch = 16,
         cex = mds.point.cex)

    # PCA for TEs ONLY
    plot(PCA.TEs$x[, 1],
         PCA.TEs$x[, 2],
         col = pca.colors,
         main = paste(treatment_label, " OE, Repeat Expression PCA", sep = ''),
         xlab = paste('PC1 (', round(100*summary.PCA.TEs$importance[,1][2],1),"%)", sep=""),
         ylab = paste('PC2 (', round(100*summary.PCA.TEs$importance[,2][2],1),"%)", sep=""),
         #xlim = c(-75, 75),
         #ylim = c(-25, 25),
         cex.main = mds.main.cex,
         cex.lab = mds.label.cex,
         cex.axis = mds.axis.cex,
         pch = 16,
         cex = mds.point.cex)

    # PCA for L1 ONLY
    plot(PCA.L1$x[, 1],
         PCA.L1$x[, 2],
         col = pca.colors,
         main = paste(treatment_label, " OE, L1 Expression PCA", sep = ''),
         xlab = paste('PC1 (', round(100*summary.PCA.L1$importance[,1][2],1),"%)", sep=""),
         ylab = paste('PC2 (', round(100*summary.PCA.L1$importance[,2][2],1),"%)", sep=""),
         #xlim = c(-25, 20),
         #ylim = c(-25, 25),
         cex.main = mds.main.cex,
         cex.lab = mds.label.cex,
         cex.axis = mds.axis.cex,
         pch = 16,
         cex = mds.point.cex)

dev.off()





  

# Return nothing
return(NULL)
  
  
  
} # END FUNCTION
