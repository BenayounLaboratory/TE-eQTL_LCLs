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