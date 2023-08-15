check_if_QTLs_replicate_function <- function(snp_info_df, genotype_BED_file_path, Expression_file_path, sig_QTL_df){
  
  # THIS FUNCTION TAKES SNPS THAT ARE SIGNIFICANT QTLS IN ONE POPULATION, AND DETERMINES THEIR SIGNIFICANCE IN ANOTHER POPULATION. IT ESSENTIALLY RUNS A LINEAR REGRESSION COMPARING GENOTYPE VS GENE EXPRESSION.
  # NOTE: THIS FUNCTION ASSUMES SNPS IN A SIGNIFICANT EQTL TABLE ARE IN A COLUMN WITH NAME 'RsID' AND THE ASSOCIATED GENE/TE IS IN A COLUMN NAMED 'gene'
  
  
  
  
  
  # Load necessary files
      
  # Load Gene/TE expression
  Gene_expression <- read.csv(Expression_file_path, header = T, row.names = 1, stringsAsFactors = F, sep = '\t')
  
  # Make R object to stream snp genotypes into memory (rows are samples and columns are snps). #NOTE: MT snps will need to be assigned a snpid, since they dont have rsid.
  binary_genotypes <- BEDMatrix(path = genotype_BED_file_path, simple_names = T)
  
      # Arrange row order (samples) to match expression df
      binary_genotypes <- binary_genotypes[colnames(Gene_expression), ]
      
      
      
  # Make dataframe to hold linear regression stats in the new population. 
  
  # Define the new df 
  new.QTLs <- sig_QTL_df[, c('RsID', 'gene')]

  # Extract custom snpIDs (these are unique; RsID may not be)
  new.QTLs$snps <- snp_info_df[sig_QTL_df$RsID, 'snpid']
    
  # Add column to hold linear regression pvalues
  new.QTLs$pvalue <- 10
  
  # Add column to hold linear regression betas
  new.QTLs$beta
  
  # Define the total number of SNV-gene interactions
  total_snp_gene_interactions <- nrow(new.QTLs)
      
  
  
  # Loop over each significant snp-gene interaction
  for (ith_interaction in 1:total_snp_gene_interactions) {
    
    
      # Define the ith gene
      ith_gene <- new.QTLs[ith_interaction, 'gene']
      
      # Define expression for the ith gene
      ith_gene_expression <- as.numeric(Gene_expression[ith_gene, ])
          
      # Define the ith snp (in RsID format, which is used in BED files)
      ith_snp <- new.QTLs[ith_interaction, 'RsID']
      
      # Define the genotypes for the ith snp
      ith_genotypes <- as.numeric(binary_genotypes[, ith_snp])
    
      # Run linear regression of expression vs snp
          
          # run regression
          lm.snp_gene_expr <- lm(ith_gene_expression ~ ith_genotypes)
      
          # Extract the linear regression pvalue
          ith_snp_gene_regression_pval <- summary(lm.snp_gene_expr)$coef["ith_genotypes","Pr(>|t|)"]
          
          # Extract the linear regression coefficient
          ith_regression_coefficient <- summary(lm.snp_gene_expr)$coef["ith_genotypes", "Estimate"]
          
          # Update the new QTL stats df with the pvalue
          new.QTLs[ith_interaction, 'pvalue'] <- ith_snp_gene_regression_pval
          
          # Update the new QTL stats df with the beta
          new.QTLs[ith_interaction, 'beta'] <- ith_regression_coefficient
    
      # Continue to the next snp-gene interaction
    
    
  } # END FOR LOOP OVER SNV-GENE INTERACTIONS
  
  
      
  # Add BH FDR values
  new.QTLs$FDR <- p.adjust(new.QTLs$pvalue, method = 'BH')
  
  # return the final QTL df
  return(new.QTLs)
  
  
} # END FUNCTION
