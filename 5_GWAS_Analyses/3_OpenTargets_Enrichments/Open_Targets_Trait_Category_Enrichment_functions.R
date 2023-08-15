get_phewas <- function(SNP_IDs_dataframe, pval_filter) {
  
  # FUNCTION INFO: RETURN CSV WITH TRAITS AND TRAIT CATEGORIES ASSOCIATED WITH EACH SNP
  # SNP_IDs_dataframe should be a dataframe containing a SNP's RsID and OpenTargets variant_id
  # pval filter is the pvalue threshold for significant associations
  
  
  
  
  # Create an empty phewas data frame
  my_phewas <- data.frame() 
  
  # Define the Open Targets variant ids to run
  if(nrow(SNP_IDs_dataframe) != 1) {
    
    my_SNP_set <- SNP_IDs_dataframe$variant_id
  }
  else {
    
    # If the dataframe has 1 row, like during the random sampling, transpose the dataframe
    my_SNP_set <- t(SNP_IDs_dataframe)
  }

  # Loop over each variant_id in my_SNP_set
  for (variant_id in my_SNP_set) {
    
      # Query Open Targets
      otg_qry$query('phewas', 'query phewas($variantInfo: String!) {
        pheWAS (variantId: $variantInfo) {
           associations {
            pval
            study {
              traitReported
              studyId
              traitCategory
            }
          } 
        }
      }')
    
      # List each variant_id as a list
      variables <- list(variantInfo = variant_id)
      
      # Flatten the JSON results of the query into something R can handle
      results <- fromJSON(otg_cli$exec(otg_qry$queries$phewas, variables, flatten = TRUE))
      
      # Create a dataframe from the results that includes the trait reported, the study ID, and the trait category
      ith_phewas <- as.data.frame(results$data$pheWAS$associations$study) 
      
      # Create a dataframe from the results that includes the pvalue
      ith_phewas_pvals <- as.data.frame(results$data$pheWAS$associations$pval) 
      
      tryCatch({
        
        # Add a new column and insert the variant ID for reference
        ith_phewas$variantInfo <- variant_id 
        
        if(nrow(SNP_IDs_dataframe) != 1) {
          # Add a new column and insert the rsid for reference
          ith_phewas$rsid <- SNP_IDs_dataframe$rsid[SNP_IDs_dataframe$variant_id == variant_id] 
        }
        
        
        # Bind the dataframe that includes trait reported, study ID, trait category, rsid, variant info, and rsid with the dataframe that consists of one column with the pvalue
        ith_phewas <- cbind(ith_phewas, ith_phewas_pvals) 
        
        # Bind the retrieved data for this one variant id with all the previous
        my_phewas <- rbind(ith_phewas, my_phewas)
      }, 
      error = function(e) {}
      )
      
  } # END FOR LOOP
  
  
  # Change the column title for the column with the p values
  colnames(my_phewas)[which(names(my_phewas) == "results$data$pheWAS$associations$pval")] <- "pval"
  
  # Filter associations on pvalue
  my_phewas <- my_phewas[which(my_phewas$pval < pval_filter), ]
  
  # Return the dataframe with phewas results for all snps tested
  return(my_phewas) 
  
  
  
} # END FUNCTION


get_p_value <- function(real_data_count, random_data_counts, num_perms) {
  
  # This function outputs a pvalue based on an ecdf distribution
  
  
  
  # Generate an ecdf function using random data
  my_ecdf <- ecdf(random_data_counts)
  
  # Calculate the raw pvalue
  p_value <- my_ecdf(real_data_count)
  
  # The raw pvalue = probability of getting value smaller than X, so this can be used in cases of depletion. 
  # To calculate enrichment on the right end, do 1-ecdf
  if(p_value > 0.5){
    p_value <- 1-p_value
  }
  if(p_value == 0) {
    p_value <- 1/num_perms
  }
  
  return(p_value)
  
  
} # END FUNCTION