generate_permutations <- function(vector_to_permute, unique_only, perm_length, max_permutations, perm_layout, perms_to_keep, output_file){
  
  # NOTE: Max permutations should be relatively high. Rows with even a single permuted value equaling the input value will be removed.
  # NOTE: Output may have the original unpermuted row/column for comparison purposes. Make sure it is read as such and not as a permutation.
  # Potential code improvement: keep either the 'colum' or 'row' loop and transpose if I want the other. This would just require specifying separate save conditions.
  # Requires the 'arrangements' library
  
  
  
  
  # Set seed for reproducibility
  set.seed(90280)
  
  # Generate a matrix of permutations
  permutation_matrix <- permutations(x = vector_to_permute, k = perm_length, nsample = max_permutations, replace = F, layout = perm_layout)
  
  
  if (perm_layout == 'column' && unique_only == T) {
    
    # Row names will reflect the non-permuted vector
    rownames(permutation_matrix) <- vector_to_permute
    
    # Loop over each sample in the unpermuted vector. IF permutations have the same sample in the same order, remove that permutation.
    for (sample_num in 1:nrow(permutation_matrix)) {
      
      # Get the ith sample name
      ith_sample_name <- rownames(permutation_matrix)[sample_num]
      
      # Keep columns (permutations) in the permutation matrix that do not have the original sample somewhere else in the same row (meaning that sample was not scrambled)
      permutation_matrix <- permutation_matrix[, which(permutation_matrix[sample_num, ] != ith_sample_name)]
      
    }
    
    # Subset the number of permutations I'm interested in keeping.
    permutation_matrix <- permutation_matrix[, 1:perms_to_keep]
    
    # Save table
    write.table(permutation_matrix, file = output_file, sep = "\t" , row.names = T, col.names = F, quote = F)
    
    # Return table (for inspection)
    return(permutation_matrix)
    
    
    
  } else if (perm_layout == 'row' && unique_only == T) {
    
    # Column names will reflect the non-permuted vector
    colnames(permutation_matrix) <- vector_to_permute
    
    # Loop over each sample in the unpermuted vector. IF permutations have the same sample in the same order, remove that permutation.
    for (sample_num in 1:ncol(permutation_matrix)) {
      
      # Get the ith sample name
      ith_sample_name <- colnames(permutation_matrix)[sample_num]
      
      # Keep rows (permutations) in the permutation matrix that do not have the original sample somewhere else in the same column (meaning that sample was not scrambled)
      permutation_matrix <- permutation_matrix[which(permutation_matrix[, sample_num] != ith_sample_name), ]
      
    }
    
    # Subset the number of permutations I'm interested in keeping.
    permutation_matrix <- permutation_matrix[1:perms_to_keep, ]
    
    # Save table
    write.table(permutation_matrix, file = output_file, sep = "\t" , row.names = F, col.names = T, quote = F)
    
    # Return table (for inspection)
    return(permutation_matrix)
    
  } else {
    
    # Return table (for inspection)
    return(permutation_matrix)
    
    # Save table
    write.table(permutation_matrix, file = output_file, sep = "\t" , row.names = F, col.names = F, quote = F)
    
  }
  
  
  
  
} #END FUNCTION
