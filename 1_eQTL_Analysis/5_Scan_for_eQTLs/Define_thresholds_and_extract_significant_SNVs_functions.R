empiricalFDR_to_pval <- function(data.file, null_data.dir, eQTL_type, arbitrary_pval, empirical_FDR_target){
  
  # eQTL_type should be 'cis' or 'trans', depending on the eQTL analysis
  # abitrary pval will be used to trim eqtl data files so that they can be loaded into memory
  # empirical FDR target is the empirical FDR we will select a corresponding pval threshold for
  # average empirical FDR at a given pvalue is defined as (number_permutation_snps/N_of_permutations)/number_real_snps
  
  # Ouput for the function below is 
  # 1) first pvalue where empirical FDR < empirical.FDR.threshold 
  # 2) the exact FDR at that pvalue 
  # 3) the average number of permutation datapoints with pvalue <= the pvalue calculated in result 1 ie false positives and 
  # 4) the number of real data points with pvalue <= the pvalue calculated in result 1. Dividing result 3 by result 4 yields the empirical FDR at that pvalue
  # Note, the function starts at 'arbitrary_pval', calculates the empirical FDR using the real and the permutation data, and consecutively goes to the next smallest pval until the target empirical FDR is reached. If it is not reached, it will spit out NULL.
 
    
    
  


  # Read eQTL raw results file
  eQTL.stats <- readRDS(data.file)
  
  # Only keep stats passing arbitrary pval threshold (to avoid loading everything into memory)
  if (eQTL_type == c('cis')) {
    
      eQTL.stats <- eQTL.stats$cis$eqtls
      eQTL.stats <- eQTL.stats[which(eQTL.stats$pvalue <= arbitrary_pval), ]
    
  } else if (eQTL_type == c('trans')) {
    
      eQTL.stats <- eQTL.stats$all$eqtls
      eQTL.stats <- eQTL.stats[which(eQTL.stats$pvalue <= arbitrary_pval), ]
    
  }
  
  
  
  # Define permutation data files
  permutations.files <- list.files(null_data.dir, "\\.rds$", recursive=FALSE, full.names=TRUE)
  
  # Initiate counter 
  permutation_counter <- 0
  
  # Aggregate stats passing the arbitrary pval threshold across ALL permutation files
  for (file in permutations.files) {
   
    # Update permutation counter
    permutation_counter <- permutation_counter + 1
    
    # Read external file 
    null.eQTL.stats <- readRDS(file)
    
    # Only keep stats (passing arbitrary pval threshold)
    if (eQTL_type == c('cis')) {
      
        null.eQTL.stats <- null.eQTL.stats$cis$eqtls
        null.eQTL.stats <- null.eQTL.stats[which(null.eQTL.stats$pvalue <= arbitrary_pval), ]
      
    } else if (eQTL_type == c('trans')) {
      
        null.eQTL.stats <- null.eQTL.stats$all$eqtls
        null.eQTL.stats <- null.eQTL.stats[which(null.eQTL.stats$pvalue <= arbitrary_pval), ]
      
    }

    # If it's the first permutation, define a new dataframe. Otherwise, append the ith permutation data to the growing dataframe.
    if (permutation_counter == 1) {
      
      all_permutation_data <- null.eQTL.stats
      
    } else {
      
      all_permutation_data <- rbind(all_permutation_data, null.eQTL.stats)
      
    }
    
     
  } # End for loop over permutation files
  
  
  
  # Order permutation data by ascending pvalue
  all_permutation_data <- all_permutation_data[order(all_permutation_data$pvalue, decreasing = FALSE), ]
  
  # Order real eQTL results by ascending pvalue
  eQTL.stats <- eQTL.stats[order(eQTL.stats$pvalue, decreasing = FALSE), ]
  
  # Start testing pvalues, going from highest to lowest (will keep the first pvalue where FDR < empirical_FDR_target)
  
  # Define the total number of pvalues to test
  total_pvalues <- nrow(eQTL.stats)
  
  # Define the total number of permutations
  N_of_permutations <- length(permutations.files)
  
  # Loop over pvalues
  for (ith_entry in (total_pvalues:1)) {
    
    # Extract the ith real pvalue
    ith_pvalue <- eQTL.stats[ith_entry, 'pvalue']
    
    # Count the number of snps with equal or smaller pvalue in the *real* data
    number_real_snps <- nrow(eQTL.stats[which(eQTL.stats$pvalue <= ith_pvalue), ])
    
    # Count the number of snps with equal or smaller pvalue in the *permutation* data
    number_permutation_snps <- nrow(all_permutation_data[which(all_permutation_data$pvalue <= ith_pvalue), ])
    
    # Obtain the average number of permutation snps (by dividing by the number of permutations)
    avg_number_permutation_snps <- number_permutation_snps/N_of_permutations
    
    # Calculate empirical FDR at the ith_pvalue
    ith_empirical_FDR <- avg_number_permutation_snps/number_real_snps
    
    # Provide an output if the target FDR is reached
    if (ith_empirical_FDR < empirical_FDR_target) {
      
      output_vector <- c(ith_pvalue, ith_empirical_FDR, avg_number_permutation_snps, number_real_snps)
      return(output_vector)
      
    }
    
  }
  
  
} # END FUNCTION

Extract_sig_snps_from_matrixEQTL <- function(output_location, eQTL_stats_df, vector_of_significance_thresholds, snp_info_df, output.file.label, save_rsid_list){
  
  # THIS FUNCTION TAKES ****PVALUES**** CORRESPONDING TO MULTIPLE SIGNIFICANT THRESHOLDS (BH FDR, EMPIRICAL FDR, GENOMEWIDE THRESHOLD, ETC.), CHOOSES THE STRICTEST ONE, EXTRACTS SNPS PASSING THAT THRESHOLD FROM A DF WITH MATRIXEQTL RESULTS, AND SAVES THEM (IF THERE IS AT LEAST ONE)
  # CAN SAVE A LIST OF UNIQUE, SIGNIFICANT RSIDS, IF DESIRED
  # IMPORTANT: IN THE CURRENT ITERATION OF THE SCRIPT, THE PVALUES CORRESPOND TO A BH OR EMPIRCAL FDR OF < 0.05. THEREFORE, INCLUDING THE PVALUES AT THE THRESHOLD WILL STILL SATISFY FDR < 0.05
  
  
  
  # Assign snpid to the rownames of the snp df for easier manipulation
  rownames(snp_info_df) <- snp_info_df$snpid
  
  # Define the strictest significance threshold
  final.significance.threshold <- min(vector_of_significance_thresholds)
  
  # Extract snps passing FDR 
  eQTL.stats.sig <- eQTL_stats_df[eQTL_stats_df$pvalue <= final.significance.threshold, ]
  
  # Add RsID to the significant snps table
  eQTL.stats.sig <- data.frame(RsID = snp_info_df[eQTL.stats.sig$snps, 'RsID'],
                                         eQTL.stats.sig)
  
  # Define unique, significant RsIDs
  unique.sig.RsIDs <- unique(eQTL.stats.sig$RsID)
  
  # Save, if there are significant snps
  if (nrow(eQTL.stats.sig) > 0) {
    
    write.table(eQTL.stats.sig, file = paste(output_location, output.file.label, Sys.Date(), '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)
    
    if (save_rsid_list == c('yes')) { # Save unique RsID list, if desired
      
       write.table(unique.sig.RsIDs, file = paste(output_location, output.file.label, 'Unique_RsIDs_', Sys.Date(), '.txt', sep =""), sep = "\t" , row.names = F, col.names = F, quote = F)
    }
   
    
  }
      
  
  
} # END FUNCTION

