eQTL_to_CMplot <- function(eQTLs_for_plot, snp_pos){
  
  # Subset the eQTL_df and define the CMplot object.
  CMplot_df <- eQTLs_for_plot[, c('snps', 'gene', 'pvalue')]
  
  # Change snp-pos rownames to make subsetting easier.
  rownames(snp_pos) <- snp_pos$snpid
  
  # Add column of unique indices (to identify unique TE-gene comparisons)
  CMplot_df$index <- rownames(CMplot_df)
  
  # Add columns with chromosome name/pos
  CMplot_df <- cbind(CMplot_df, snp_pos[CMplot_df$snps, c('chr', 'pos')])
  
  # Rearrange columns
  CMplot_df <- CMplot_df[, c('index', 'chr', 'pos', 'pvalue', 'snps', 'gene')]
  
  # Return the final object to plot
  return(CMplot_df)
  
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

eQTL_linear_regression_gene_TE <- function(linked_cis_trans_QTL_df, expression_df){
  
  # THIS FUNCTION TAKES A DF WITH GENES AND TES LINKED BY A COMMON SNP, RUNS A LINEAR REGRESSION ON THE EXPRESSION OF THE GENE AND THE TE, AND ADDS THE PVALUE FOR THE GENE-TE RELATIONSHIP TO THE DF (TO ASSESS SIGNIFICANCE OF THE RELATIONSHIP)
  
  
  

  
  # Define the total number of potential snp-gene-TE relationships
  total_interactions <- nrow(linked_cis_trans_QTL_df)
  
  # Add column for regression pvalue
  linked_cis_trans_QTL_df$Gene_TE_pval <- 10
  
  # Add column for regression FDR
  linked_cis_trans_QTL_df$Gene_TE_FDR <- 10
  
  # Add column for regression coefficient
  linked_cis_trans_QTL_df$Gene_TE_beta <- 10
  
  # loop over each interaction
  for (ith_relationship in 1:total_interactions) {
    
      # if snp was not linked to a gene, skip to the next loop iteration
      if (is.na(linked_cis_trans_QTL_df[ith_relationship, 'gene'])) {
            next
      }
    
      # Extract labels for the ith gene-TE snp
        
          # Extract the ith relationship gene label
          gene_i <- linked_cis_trans_QTL_df[ith_relationship, 'gene']
          
          # Extract the ith relationship TE label
          TE_i <- linked_cis_trans_QTL_df[ith_relationship, 'TE']
            

      # Extract expression measures
          
          # gene expression
          gene_i_expression <- as.numeric(expression_df[gene_i,])
          
          # TE expression
          TE_i_expression <- as.numeric(expression_df[TE_i, ])
    
          
      # linear regression to compare TE expression vs gene expression
                
          # run regression
          lm.gene_TE_expr <- lm(TE_i_expression ~ gene_i_expression)

          # Extract the linear regression pvalue
          ith_regression_pval <- summary(lm.gene_TE_expr)$coef["gene_i_expression","Pr(>|t|)"]
          
          # Extract the linear regression coefficient
          ith_regression_coefficient <- summary(lm.gene_TE_expr)$coef["gene_i_expression", "Estimate"]
          
          
      # Update the original df with regression pval
      linked_cis_trans_QTL_df[ith_relationship, 'Gene_TE_pval'] <- ith_regression_pval
      
      # Update the original df with regression coefficient
      linked_cis_trans_QTL_df[ith_relationship, 'Gene_TE_beta'] <- ith_regression_coefficient
      
    
      # Continue to next loop iteration
    
    
  } # End for loop
  
  
  
  
  # For any unmapped relationships, add NA for the pvalue, FDR, and coefficient
  linked_cis_trans_QTL_df[which(linked_cis_trans_QTL_df$Gene_TE_pval == 10), c('Gene_TE_pval', 'Gene_TE_FDR', 'Gene_TE_beta')] <- NA

  # Add linear regression BH FDR
  linked_cis_trans_QTL_df$Gene_TE_FDR <- p.adjust(linked_cis_trans_QTL_df$Gene_TE_pval, method = 'BH')
  

  # Return the updated df
  return(linked_cis_trans_QTL_df)
  
  
  
} # END FUNCTION

link_cis_and_trans_eQTLs <- function(cis.sig, trans.sig, gene_TE_INT_expression){
  
  # THIS FUNCTION: TAKES A TABLE WITH TRANS-EQTLS AND ONE WITH CIS-EQTLS, AND MAKES A ROW FOR EACH SNP-TE-GENE TRIO INTERACTION. THEN, USING AN INPUT EXPRESSION DATAFRAME, A LINEAR REGRESSION IS CARRIED OUT TO COMPARE GENE AND TE EXPRESSION
  # USES ENSEMBL_TO_SYMBOL FUNCTION
  # USES eQTL_linear_regression_gene_TE FUNCTION
  
  
  

  # Define counter to keep track of TE snps 
  counter <- 0

  # Order cis snps by pvalue (smallest to largest)
  cis.sig <- cis.sig[order(cis.sig$pvalue, decreasing = FALSE), ]
  
  # Run loop over each snp in the significant TE trans-eQTL table, and combine trans-eQTL and cis-eQTL stats
  for (TE_snp_index in 1:nrow(trans.sig)) {
    
    counter <- counter + 1
    
    # Define the ith TE snp
    ith_TE_snp <- trans.sig[TE_snp_index, 'snps']
    
    # Define the gene stats associated with that snp
    ith_TE_snp_gene_stats <- cis.sig[cis.sig$snps == ith_TE_snp,]
    
    # if TE snp is not a gene SNP, add placeholder stats. otherwise, add the gene stats
    if (nrow(ith_TE_snp_gene_stats) == 0) {
      
      # Make temporary, dummy df to hold gene AND TE stats for the ith SNP 
      TE_Gene_stats <- as.data.frame(matrix(10, nrow = 1, ncol = 10))
      colnames(TE_Gene_stats) <- c('RsID', 'snps', 'gene', 'gene_pvalue', 'gene_FDR', 'gene_beta', 'TE', 'TE_pvalue', 'TE_FDR', 'TE_beta')
      
      # Assign TE stats to temporary df
      TE_Gene_stats[, c('RsID', 'snps', 'TE', 'TE_pvalue', 'TE_FDR', 'TE_beta')] <- trans.sig[TE_snp_index, c('RsID', 'snps' ,'gene', 'pvalue', 'FDR', 'beta')]
      

    } else {
      
      # Make temporary, empty df to hold gene AND TE stats for the ith SNP 
      TE_Gene_stats <- as.data.frame(matrix(10, nrow = nrow(ith_TE_snp_gene_stats), ncol = 10))
      colnames(TE_Gene_stats) <- c('RsID', 'snps', 'gene', 'gene_pvalue', 'gene_FDR', 'gene_beta', 'TE', 'TE_pvalue', 'TE_FDR', 'TE_beta')
      
      # Assign gene stats to temporary df
      TE_Gene_stats[, c('gene', 'gene_pvalue', 'gene_FDR', 'gene_beta')] <- ith_TE_snp_gene_stats[, c('gene', 'pvalue', 'FDR', 'beta')]
      
      # Assign TE stats to temporary df
      TE_Gene_stats[, c('RsID', 'snps', 'TE', 'TE_pvalue', 'TE_FDR', 'TE_beta')] <- trans.sig[TE_snp_index, c('RsID', 'snps', 'gene', 'pvalue', 'FDR', 'beta')]
    
    
    }
    
    
    
    # If this is the first snp, make a df to hold the results. Otherwise, add results to the existing dataframe
    if (counter == 1) {
      ALL_TE_Gene_stats <- TE_Gene_stats
    } else if (counter > 1) {
      ALL_TE_Gene_stats <- rbind(ALL_TE_Gene_stats, TE_Gene_stats)
    }
    
    
    
  } # End for loop over TE snps
  
  
  
  
  
  # ADD GENE SYMBOL INFO
  
  # Define indices for *actually* linked Genes (ie that do not have the placeholder)
  actually_linked_gene_indices <- which(ALL_TE_Gene_stats$gene != 10)
  
  # Define EnsemblIDs for *actually* linked Genes (ie that do not have the placeholder)
  actually_linked_genes <- ALL_TE_Gene_stats[actually_linked_gene_indices, 'gene']
  
  # Add column to hold gene symbols or names
  ALL_TE_Gene_stats$symbol <- 10
  
  # if there were no linked genes, do nothing. otherwise, obtain gene names
  if (length(actually_linked_gene_indices) == 0) {
    
      # Do nothing
    
  } else {
      
      # Run function that maps EnsemblID to symbol
      Mapped_gene_info <- Ensembl_to_symbol(organism = 'hs', Ensembl_genes = actually_linked_genes)

      # Extract hgnc symbols for genes of interest
      ALL_TE_Gene_stats[actually_linked_gene_indices, 'symbol'] <- Mapped_gene_info[actually_linked_genes, 'hgnc_symbol'] 
    
  }
  
        
  
  
  
  # CLEANUP THE DF
  
  # Change all unmapped gene info to NA
  ALL_TE_Gene_stats[which(ALL_TE_Gene_stats$gene == 10), c('gene', 'gene_pvalue', 'gene_FDR', 'gene_beta', 'symbol')] <- NA
        
  # Order results by TE trans-eQTL pvalue 
  ALL_TE_Gene_stats <- ALL_TE_Gene_stats[order(ALL_TE_Gene_stats$TE_pvalue, decreasing = FALSE), ]
  
  # Reorder columns
  ALL_TE_Gene_stats <- ALL_TE_Gene_stats[, c('snps', 'RsID', 'gene', 'symbol', 'gene_pvalue', 'gene_FDR', 'gene_beta', 'TE', 'TE_pvalue', 'TE_FDR', 'TE_beta')]
  

  
  
  
  # RUN LINEAR REGRESSION (COMPARING GENE AND TE EXPRESSION)
  ALL_TE_Gene_stats <- eQTL_linear_regression_gene_TE(linked_cis_trans_QTL_df = ALL_TE_Gene_stats, expression_df = gene_TE_INT_expression)
  
  
  
  
  
  # Return the combined dataframe 
  return(ALL_TE_Gene_stats)
  
  
  
} # END FUNCTION

eQTL_to_barplots_scatterplot <- function(output.dir, output.filename, plot_rows, plot_columns, indices_to_plot, combined_cis_trans_QTL_df, all_snp_positions, gene_expression_df, streamed_genotypes_matrix){
  
# THIS FUNCTION WILL TAKE A TABLE WITH LINKED SNP-GENE-TE RELATIONSHIPS AND PRODUCE 3 PLOTS (SNP-GENE BARPLOT, SNP-TE BARPLOT, GENE-TE SCATTERPLOT), FOR AT MOST 3 SNPS
# REQUIRES THE BEDMATRIX library




# Save as PDF
pdf(paste(output.dir, output.filename, Sys.Date(), ".pdf", sep=""), width = 8, height = 8)
par(mfrow = c(plot_rows, plot_columns), pty ="s") # Produce plot with specified rows and columns

    # Loop over the indices for the snp-gene-TE pairs you'd like to plot
    for (ith_gene_TE_pair in indices_to_plot) {
      
      
      # EXTRACT GENOTYPES, EXPRESSION VALUES, AND LABELS
      
      # Extract labels for the ith gene-TE snp
      
          # Extract the ith gene label
          gene_i <- combined_cis_trans_QTL_df[ith_gene_TE_pair, 'gene']
          
          # Extract the ith gene symbol
          gene_i_symbol <- combined_cis_trans_QTL_df[ith_gene_TE_pair, 'symbol']
          
          # Extract the ith TE label
          TE_i <- combined_cis_trans_QTL_df[ith_gene_TE_pair, 'TE']
          
          # Extract the ith snpid
          snp_i <- combined_cis_trans_QTL_df[ith_gene_TE_pair, 'snps']
          
          # Extract the ith RsID
          snp_i_RSid <- combined_cis_trans_QTL_df[ith_gene_TE_pair, 'RsID']
      
      
      # Extract expression values
          
          # gene expression
          gene_i_expression <- gene_expression_df[gene_i,]
          
          # TE expression
          TE_i_expression <- gene_expression_df[TE_i, ]
          
      # Extract genotypes for the ith snp
          
          # In 0/1/2 format
          snp_i_genotypes <- streamed_genotypes_matrix[colnames(gene_expression_df), snp_i_RSid]

          # assign colors to each genotype (used to specify the point of each color in the gene-TE scatterplot)
          snp_i_genotypes.colors <- snp_i_genotypes
          snp_i_genotypes.colors[snp_i_genotypes.colors == 0] <- 'coral'
          snp_i_genotypes.colors[snp_i_genotypes.colors == 1] <- 'firebrick2'
          snp_i_genotypes.colors[snp_i_genotypes.colors == 2] <- 'firebrick4'
          
          # define the three genotypes by A/C/T/G allele
          snp_i.alleles.0 <- paste(all_snp_positions[snp_i, 'REF'], all_snp_positions[snp_i, 'REF'], sep = '')
          snp_i.alleles.1 <- paste(all_snp_positions[snp_i, 'REF'], all_snp_positions[snp_i, 'ALT'], sep = '')
          snp_i.alleles.2 <- paste(all_snp_positions[snp_i, 'ALT'], all_snp_positions[snp_i, 'ALT'], sep = '')
          
                # collect the 3 genotypes by allele (used for more specific x-axis labeling instead of just 0/1/2)
                snp_i.alleles <- c(snp_i.alleles.0, snp_i.alleles.1, snp_i.alleles.2)
                
                # assign colors to the three genotypes allele pairs (used to specify colors in the legends)
                snp_i.alleles.colors <- c('coral', 'firebrick2', 'firebrick4')
                
                
      
      # RUN LINEAR REGRESSIONS
                
      # make df with all ith snp-gene-info
      all_data_ith_comparison <- data.frame(genotype = snp_i_genotypes,
                                            gene_expr = as.numeric(gene_i_expression),
                                            TE_expr = as.numeric(TE_i_expression)
                                            )
      
      # copy the df !!!!!!!!!!!!!!!!!!!! BUT with a +1 offset for the genotypes 
      # this is necessary since since boxplots are generated for 1st, 2nd...nth categories. this will be used to shift the corresponding snp-gene/snp-TE linear regression intercept so it can be plotted appropriately on the boxplots.
      all_data_ith_comparison.Genotype_OFFSET <- all_data_ith_comparison
      all_data_ith_comparison.Genotype_OFFSET$genotype <- all_data_ith_comparison.Genotype_OFFSET$genotype + 1
      
      # run linear regression to compare snps vs Gene expression
          
          # run regression
          lm.snp_gene_expr <- lm(gene_expr ~ genotype, data = all_data_ith_comparison.Genotype_OFFSET)
      
          # Extract the linear regression pvalue
          ith_snp_gene_regression_pval <- summary(lm.snp_gene_expr)$coef["genotype","Pr(>|t|)"]
          
          # Only keep two pval decimals (pval is included in the scatterplot legend)
          ith_snp_gene_regression_pval <- formatC(ith_snp_gene_regression_pval, format = 'e', digits = 2)
          
      # linear regression to compare snps vs TE expression
          
          # run regression
          lm.snp_TE_expr <- lm(TE_expr ~ genotype, data = all_data_ith_comparison.Genotype_OFFSET)
      
          # Extract the linear regression pvalue
          ith_snp_TE_regression_pval <- summary(lm.snp_TE_expr)$coef["genotype","Pr(>|t|)"]
          
          # Only keep two pval decimals (pval is included in the scatterplot legend)
          ith_snp_TE_regression_pval <- formatC(ith_snp_TE_regression_pval, format = 'e', digits = 2)
          
      # linear regression to compare TE expression vs gene expression
          
          # run regression
          lm.gene_TE_expr <- lm(TE_expr ~ gene_expr, data = all_data_ith_comparison)
      
          # Extract the linear regression pvalue
          ith_regression_pval <- summary(lm.gene_TE_expr)$coef["gene_expr","Pr(>|t|)"]
          
          # Only keep two pval decimals (pval is included in the scatterplot legend)
          ith_regression_pval <- formatC(ith_regression_pval, format = 'e', digits = 2)
  
      
          
      # GENERATE THE PLOTS
              
      # TE plot
          
          # make boxplots for genotype segregated TE expression
          boxplot(TE_expr ~ genotype,
                  data = all_data_ith_comparison,
                  ylim = c(-4, 4),
                  outline = F,
                  #main = paste(TE_i, 'vs', snp_i_RSid), 
                  xlab = 'Genotype', 
                  ylab = paste('INT(', TE_i, ')', sep = ''),
                  xaxt = 'n')

          # change title so it is not bold
          # title(main = paste(TE_i, 'vs', snp_i_RSid), font.main = 1, cex.main = 1)

          # overlay the points
          stripchart(TE_expr ~ genotype, data = all_data_ith_comparison, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 0.5, col = snp_i.alleles.colors) # Change axis to fit all points

          # add linear regression line
          abline(lm.snp_TE_expr, col = 'black', lty = 1)

          # legend with regression pval
          legend("top", c(ith_snp_TE_regression_pval), col = c('white'), pch = 16, cex = 1, pt.cex = 0, box.lty = 1, box.lwd = 0, box.col = 'white', bg = 'transparent')

          # change 0/1/2 genotype axis labels to allele pair labels
          axis(1, at = 1:3, snp_i.alleles)
      
          
      
      # Gene plot
          
          # make boxplots for genotype segregated gene expression
          boxplot(gene_expr ~ genotype, 
                  data = all_data_ith_comparison,
                  ylim = c(-4, 4),
                  outline = F, 
                  #main = paste(gene_i_symbol, 'vs', snp_i_RSid), 
                  xlab = 'Genotype', 
                  ylab = paste('INT(', gene_i_symbol, ')', sep = ''),
                  xaxt = 'n')
          
          # change title so it is not bold
          #title(main = paste(gene_i_symbol, 'vs', snp_i_RSid), font.main = 1, cex.main = 1)

          # overlay the points
          stripchart(gene_expr ~ genotype, data = all_data_ith_comparison, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 0.5, col = snp_i.alleles.colors) # Change axis to fit all points
          
          # add linear regression line
          abline(lm.snp_gene_expr, col = 'black', lty = 1)

          # legend with regression pval
          legend("top", c(ith_snp_gene_regression_pval), col = c('white'), pch = 16, cex = 1, pt.cex = 0, box.lty = 1, box.lwd = 0, box.col = 'white', bg = 'transparent')
          
          # change 0/1/2 genotype axis labels to allele pair labels
          axis(1, at = 1:3, snp_i.alleles)        
          
          
          
      # Gene/TE Expression Scatter Plot
          
          # scatterplot gene  and TE expression
          plot(x = all_data_ith_comparison$gene_expr, 
               y = all_data_ith_comparison$TE_expr, 
               xlim = c(-4, 4),
               ylim = c(-4, 4),
               xlab = paste('INT(', gene_i_symbol, ')', sep = ''), 
               ylab = paste('INT(', TE_i, ')', sep = ''), 
               #main = paste(gene_i_symbol, 'vs', TE_i), 
               col = snp_i_genotypes.colors, 
               pch = 19, 
               cex = 0.50)
          
          # change title so it is not bold
          #title(main = paste(gene_i_symbol, 'vs', TE_i), font.main = 1, cex.main = 1)
          
          # Add legend with lienar regression pvalue and allelic genotope coloring
          legend("top", c(ith_regression_pval), col = c('white'), pch = 16, cex = 1, pt.cex = 0, box.lty = 1, box.lwd = 0, box.col = 'white', bg = 'transparent')
        
          # add linear regression line
          abline(lm.gene_TE_expr, col = 'black', lty = 1)
          
          # addx=0 and y=0 ticks
          abline(h = 0, lty = 2, col = 'grey60')
          abline(v = 0, lty = 2, col = 'grey60')
          

    } # End for loop over snp-gene-TE indices


dev.off()



return()
  
  
} # END FUNCTION

