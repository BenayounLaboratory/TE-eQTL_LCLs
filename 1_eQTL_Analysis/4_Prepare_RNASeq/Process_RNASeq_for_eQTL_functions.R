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

prefilter_TE_aggregation <- function(gene_TE_counts, aggregate_all, intergenic_nearby_repeats, intergenic_distal_repeats, intronic_repeats, exonic_repeats) {
  
  # FUNCTION INFO: THIS FUNCTION AGGREGATES TE COUNTS BY the specified TE group (subfamily, family, class).
  # Aggregation with this function is prior to the filtering step, since filtering removes too many loci


  
  
  # Assign name to column with gene/TE names.
  names(gene_TE_counts)[1] <- 'TEs'  
  
  # Subset genes/EBV
  gene_counts <- gene_TE_counts[grepl('ENSG', gene_TE_counts[,1]) | grepl('EBV', gene_TE_counts[,1]), ]
  
  # Subset TEs (shouldn't have ensembl prefix or be EBV expression)
  TE_counts <- gene_TE_counts[!grepl('ENSG', gene_TE_counts[,1]) | grepl('EBV', gene_TE_counts[,1]), ]
  
      # Print out how many TEs remain
      message(paste(nrow(TE_counts), 'TE loci with at least 1 count in 1 sample'))
  
  # Split TE column. Split columns are placed at end of the matrix.
  TE_counts <- as.data.frame(splitstackshape::cSplit(TE_counts, splitCols = "TEs", sep = ":", direction = "wide", fixed = TRUE, drop = TRUE, stripWhite = FALSE, makeEqual = TRUE))
  
      # Re-add the TE columns to the start of the df and delete the TE columns at the end of the df
      TE_counts <- cbind(TE_counts[, c('TEs_1', 'TEs_2', 'TEs_3', 'TEs_4')], TE_counts[, 1:(ncol(TE_counts)-4)])
    
      # Re-label TE columns
      names(TE_counts)[1:4] <- c('Locus', 'Subfamily', 'Family', 'Class')
    
      # Define label indices (columns that have label)
      label_cols <- c(1:4)
      
      # Assign locus name to the rownames
      rownames(TE_counts) <- TE_counts$Locus
      
      
  # Subset loci of interest and aggregate by subfamily
  if (aggregate_all == TRUE) { # aggregate all loci without filtering
    
      # Generate new names in format Subfamily:Family:Class
      final_names <- paste(TE_counts$Subfamily, TE_counts$Family, TE_counts$Class, sep = ':')
      
      # Aggregate counts by subfamily
      TE_counts.family.sum <- aggregate(x = TE_counts[, -label_cols], by = list(final_names), FUN = 'sum')
      
      # Name TE column; rbind doesn't work without this.
      names(TE_counts.family.sum)[1] <- 'TEs' 
      
      # Recombine gene and TE counts
      gene_TE_counts.sum <- rbind(gene_counts, TE_counts.family.sum)
    
  } else {
    
    # Group expressed repeats by genic location
    TE_counts.intergenic.nearby <- TE_counts[which(rownames(TE_counts) %in% rownames(intergenic_nearby_repeats)), ]
    TE_counts.intergenic.distal <- TE_counts[which(rownames(TE_counts) %in% rownames(intergenic_distal_repeats)), ]
    TE_counts.intronic <- TE_counts[which(rownames(TE_counts) %in% rownames(intronic_repeats)), ]
    TE_counts.exonic <- TE_counts[which(rownames(TE_counts) %in% rownames(exonic_repeats)), ]
    
    # Generate new names in format group:Subfamily:Family:Class
    intergenic.nearby.names <- paste('intergenic_near', TE_counts.intergenic.nearby$Subfamily, TE_counts.intergenic.nearby$Family, TE_counts.intergenic.nearby$Class, sep = ':')
    intergenic.distal.names <- paste('intergenic_distal', TE_counts.intergenic.distal$Subfamily, TE_counts.intergenic.distal$Family, TE_counts.intergenic.distal$Class, sep = ':')
    intronic.names <- paste('intronic', TE_counts.intronic$Subfamily, TE_counts.intronic$Family, TE_counts.intronic$Class, sep = ':')
    exonic.names <- paste('exonic', TE_counts.exonic$Subfamily, TE_counts.exonic$Family, TE_counts.exonic$Class, sep = ':')
    
    # Aggregate counts by subfamily
    TE_counts.intergenic.nearby <- aggregate(x = TE_counts.intergenic.nearby[, -label_cols], by = list(intergenic.nearby.names), FUN = 'sum')
    TE_counts.intergenic.distal <- aggregate(x = TE_counts.intergenic.distal[, -label_cols], by = list(intergenic.distal.names), FUN = 'sum')
    TE_counts.intronic <- aggregate(x = TE_counts.intronic[, -label_cols], by = list(intronic.names), FUN = 'sum')
    TE_counts.exonic <- aggregate(x = TE_counts.exonic[, -label_cols], by = list(exonic.names), FUN = 'sum')
    
    # Name the TE column; rbind doesn't work without this.
    names(TE_counts.intergenic.nearby)[1] <- 'TEs' 
    names(TE_counts.intergenic.distal)[1] <- 'TEs' 
    names(TE_counts.intronic)[1] <- 'TEs' 
    names(TE_counts.exonic)[1] <- 'TEs' 
    
    # Recombine gene and TE counts
    gene_TE_counts.sum <- rbind(gene_counts, TE_counts.intergenic.nearby, TE_counts.intergenic.distal, TE_counts.intronic, TE_counts.exonic)
  
  }


      
  
  # Return the updated counts
  return(gene_TE_counts.sum)

      
  
} # END FUNCTION

cleanup_and_filter_counts <- function(count_data, filter, min_counts_per_sample, fraction_of_samples) {
  
  # FUNCTION INFO:
  # CLEANUP: REMOVE GENE/TRANSCRIPT VERSION AND COMBINE COUNTS FROM THE SAME GENE
  # FILTER: LOWLY EXPRESSED GENES
  
  
  
  # Assign name to the column with gene/TE names. Name needs to specifically be 'gene.TE' since that is referenced downstream.
  names(count_data)[1] <- 'gene.TE'    

  
  # Remove transcript or ENSEMBL version info. 
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

Inverse_normal_transform <- function(input.df){
  # This function generates gene-wise INT expression values 
  
  # This function relies on the RankNorm function from RNOmni package
  
  
  # Matrix to hold INT values
  input.df.INT <- input.df * 0

  # Run loop that calculates INT on each gene
  for (ith_gene in 1:nrow(input.df)) {
    
    input.df.INT[ith_gene, ] <- RankNorm(as.numeric(input.df[ith_gene, ]))
    
  }
  
  
  # Output
  return(input.df.INT)
  
  
} # END FUNCTION

Split_and_save_expression <- function(output.dir, output_file_prefix, expression_df, gene_prefix, TE_groups){
  # FUNCTION INFO: This function splits an expression dataframe into All Genes/TEs, Genes only, L1s only, and Alus only expression files.
  
  
  
  
  # Check if output directory exists. If not, make it. 
  if ( dir.exists(output.dir) == FALSE) {
    dir.create(output.dir, recursive = TRUE)
  }
  
  
  # Define TE families/subfamilies
  if ( TE_groups == 'Family' ) {
    L1_prefix <- 'L1'
    Alu_prefix <- 'Alu'
  } else if ( TE_groups == 'Subfamily' ) {
    L1_prefix <- ':L1:'
    Alu_prefix <- ':Alu:'
  }
  
  
  # All genes/TEs
  write.table(expression_df, file = paste(output.dir, output_file_prefix, "All_Genes_TEs", '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)

  # Only genes
  expression_Genes <- expression_df[grepl(gene_prefix, rownames(expression_df)), , drop = FALSE]
  write.table(expression_Genes, file = paste(output.dir, output_file_prefix, "Only_Genes", '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)
    
  # All L1s
  expression_L1 <- expression_df[grepl(L1_prefix, rownames(expression_df)), , drop = FALSE]
  write.table(expression_L1, file = paste(output.dir, output_file_prefix, "Only_L1", '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)
   
  # All Alu
  expression_Alu <- expression_df[grepl(Alu_prefix, rownames(expression_df)), , drop = FALSE]
  write.table(expression_Alu, file = paste(output.dir, output_file_prefix, "Only_Alu", '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)

  
  
} # END FUNCTION

Split_and_save_expression_TElocal_quant <- function(output.dir, output_file_prefix, expression_df, gene_prefix, TE_groups){
  # FUNCTION INFO: This function splits an expression dataframe into All Genes/TEs, Genes only, L1s only, and Alus only expression files.
  
  
  
  
  # Check if output directory exists. If not, make it. 
  if ( dir.exists(output.dir) == FALSE) {
    dir.create(output.dir, recursive = TRUE)
  }
  
  
  # Define TE families/subfamilies
  if ( TE_groups == 'Family' ) {
    
    L1_prefix <- 'L1'
    
  } else if ( TE_groups == 'Subfamily' ) {
    
    L1_prefix <- ':L1:'
    
  }
  
  
  # All genes/TEs
  write.table(expression_df, file = paste(output.dir, output_file_prefix, "All_Genes_TEs", '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)

  # Only genes
  expression_Genes <- expression_df[grepl(gene_prefix, rownames(expression_df)), , drop = FALSE]
  write.table(expression_Genes, file = paste(output.dir, output_file_prefix, "Only_Genes", '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)
    
  # All L1s
  
      # Subset all L1 expression
      expression_L1 <- expression_df[grepl(L1_prefix, rownames(expression_df)), , drop = FALSE]
      write.table(expression_L1, file = paste(output.dir, output_file_prefix, "Only_L1", '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)
      
      # Subset all L1 expression (intergenic, near)
      expression_L1_intergenic_near <- expression_L1[grepl('intergenic_near', rownames(expression_L1)), , drop = FALSE]
      write.table(expression_L1_intergenic_near, file = paste(output.dir, output_file_prefix, "Only_L1_intergenic_near", '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)
      
      # Subset all L1 expression (intergenic, distal)
      expression_L1_intergenic_distal <- expression_L1[grepl('intergenic_distal', rownames(expression_L1)), , drop = FALSE]
      write.table(expression_L1_intergenic_distal, file = paste(output.dir, output_file_prefix, "Only_L1_intergenic_distal", '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)
      
      # Subset all L1 expression (intronic)
      expression_L1_intronic <- expression_L1[grepl('intronic', rownames(expression_L1)), , drop = FALSE]
      write.table(expression_L1_intronic, file = paste(output.dir, output_file_prefix, "Only_L1_intronic", '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)
      
      # Subset all L1 expression (exonic)
      expression_L1_exonic <- expression_L1[grepl('exonic', rownames(expression_L1)), , drop = FALSE]
      write.table(expression_L1_exonic, file = paste(output.dir, output_file_prefix, "Only_L1_exonic", '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)

  
  
} # END FUNCTION
