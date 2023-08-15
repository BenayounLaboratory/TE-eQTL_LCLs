load_eQTL_genotypes <- function(expression_file, SNP_file){
  
  # THIS FUNCTION LOADS THE GENOTYPES IN A FORMAT REQUIRED BY MATRIXEQTL
  # INPUTS ARE AN INT GENE EXPRESSION FILE AND THE 0/1/2 GENOTYPE MATRIX
  
  
  
  
  
  # Load gene expression data
  eqtl_phenotype <- read.csv(expression_file, header = T, row.names = 1, stringsAsFactors = F, sep = '\t')  
  
  # Use the gene expression sample order across covariates table, genotypes, and gene expression.
  eqtl_sample_order <- colnames(eqtl_phenotype) 
  
  # Remove the gene expression data
  rm(eqtl_phenotype)
  
  
  # Load genotype data in matrixeQTL format
  snps = SlicedData$new()
  snps$fileDelimiter = "\t"      # the TAB character
  snps$fileOmitCharacters = 'NA'   # denote missing values
  snps$fileSkipRows = 1          # one row of column labels
  snps$fileSkipColumns = 1       # one column of row labels
  snps$fileSliceSize = 20000      # read file in slices of X rows
  snps$LoadFile(SNP_file)
  
    # Extract and order the sample indices to match the gene expression order
    genotype_indices <- match(eqtl_sample_order, snps$columnNames)
  
  # Reorder/subset the samples in the desired order
  snps$ColumnSubsample(genotype_indices)
  
  
  
  # Return Matrix_eQTL snps object
  return(snps)

  
} # END FUNC

Single_Matrix_eQTL <- function(is.perm, output_file_cis, snps_pos, gene_pos, pvOutputThreshold_cis, cis_Dist, 
                               output_file_tra, pvOutputThreshold_tra, expression_file, gene, SNP_file, snps, covariates_file, 
                               R_obj, output.dir, R_object_filename){
  
  # Function info: This function carries out one eQTL analysis using Matrix-eQTL. This function is called by another function, Permuted_eQTLs, to run multiple eQTL analyses on permuted data.
    # SNPs can be loaded inside this function each time it's called, or they can be loaded once independently of this function (useful if the same genotypes will be used in multiple eQTL analyses)
    # cisDist is distance for local gene-SNP pairs
    # Set pv_cis=0 if you only want trans-eQTL or vice versa.
  
  
  
  
  
  # Check if output directory exists. If not, make it. 
  if ( dir.exists(output.dir) == FALSE) {
    dir.create(output.dir, recursive = TRUE)
  }
  
  
 
  
  
  # Linear model to use: modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  useModel = modelLINEAR

  # Error covariance matrix. Set to numeric() for identity.
  errorCovariance = numeric()
  
  # Load snp position file
  snps_pos <- read.table(snps_pos, header = TRUE, stringsAsFactors = FALSE)
  snps_pos <- snps_pos[, 1:3] # Input needs to have only 3 columns
  
  # Load gene position file
  gene_pos <- read.table(gene_pos, header = TRUE, stringsAsFactors = FALSE)
    
  # Load gene expression data
  eqtl_phenotype <- read.csv(expression_file, header = T, row.names = 1, stringsAsFactors = F, sep = '\t')  
  
  # Use the gene expression sample order across covariates table, genotypes, and gene expression.
  eqtl_sample_order <- colnames(eqtl_phenotype) 
  
  # Remove gene expression data
  rm(eqtl_phenotype)
  
  
  ## Load genotype data (if not already loaded outside of this function) in matrixeqtl format
  
    if (is.null(SNP_file)) {
      # Do nothing. "snps" object should be defined in input
    } else {
      # Run function to load genotypes
      snps <- load_eQTL_genotypes(expression_file = expression_file, SNP_file = SNP_file)
    }
    
  
  ## Load gene expression data in matrixeqtl format
    
    if (is.perm == TRUE) {
      # Do nothing if a permutation analysis is being run. "gene" will be specified as part of the permutation analysis function.
    } else {
      gene = SlicedData$new()
      gene$fileDelimiter = "\t"      # the TAB character
      gene$fileOmitCharacters = "NA" # denote missing values;
      gene$fileSkipRows = 1          # one row of column labels
      gene$fileSkipColumns = 1       # one column of row labels
      gene$fileSliceSize = 2000      # read file in slices of 2,000 rows
      gene$LoadFile(expression_file)
      
      # Define samples indices in the order that matches sample order in the gene expression file (in this case, the indices should already be in order)
      expression_indices <- match(eqtl_sample_order, gene$columnNames)
      
      # Reorder/subset the samples in the desired order
      gene$ColumnSubsample(expression_indices)
    }
    

    
  ## Load covariates in matrixeqtl format
    
    if (is.null(covariates_file)) {
      
      cvrt = SlicedData$new()
      
    } else {
      
      cvrt = SlicedData$new()
      cvrt$fileDelimiter = "\t"      # the TAB character
      cvrt$fileOmitCharacters = "NA" # denote missing values;
      cvrt$fileSkipRows = 1          # one row of column labels
      cvrt$fileSkipColumns = 1       # one column of row labels
      cvrt$LoadFile(covariates_file)
      
      # Define samples indices in the order that matches sample order in the gene expression file
      cov_indices <- match(eqtl_sample_order, cvrt$columnNames)
      
      # Reorder/subset the samples in the desired order
      cvrt$ColumnSubsample(cov_indices)
      
    }
  
  
  ## Run eQTL analysis using Matrixeqtl
  My_eQTL_results = Matrix_eQTL_main(output_file_name = output_file_tra,
                                     pvOutputThreshold = pvOutputThreshold_tra,
                                     output_file_name.cis = output_file_cis,
                                     pvOutputThreshold.cis = pvOutputThreshold_cis,
                                     cisDist = cis_Dist,
                                     snpspos = snps_pos,
                                     genepos = gene_pos,
                                     useModel = useModel,
                                     snps = snps,
                                     gene = gene,
                                     cvrt = cvrt,
                                     errorCovariance = errorCovariance,
                                     pvalue.hist = "qqplot",
                                     min.pv.by.genesnp = FALSE,
                                     noFDRsaveMemory = FALSE,
                                     verbose = TRUE)
    
  ## Extract statistic tables and save
  if ( pvOutputThreshold_tra == 0 ) {
    # this is for cis-eQTLs
    
    # Save cis-eQTL object
    saveRDS(object = My_eQTL_results, file = paste(output.dir, R_object_filename, ".rds", sep = ""), compress = TRUE)

          
  } else {
    # this is for trans-eQTLs
    
    # Save trans-eQTL object
    saveRDS(object = My_eQTL_results, file = paste(output.dir, R_object_filename, ".rds", sep = ""), compress = TRUE)

  }
    
  
  
  
  # Function output
  if (R_obj == TRUE) {
      return(My_eQTL_results)
  } else {
      rm(My_eQTL_results)
  }

  
} # END FUNCTION

Permuted_eQTLs <- function(permutations_file, specific_permutations, 
                           output_file_cis, snps_pos, gene_pos, pvOutputThreshold_cis, cis_Dist, 
                           output_file_tra, pvOutputThreshold_tra, expression_file, snps, covariates_file, 
                           output.dir, R_object_filename){
  
  # FUNCTION INFO: Runs Matrix_eQTL with permuted gene expression values and outputs 1) MatrixeQTL R object for each permutation, 2) vector with lowest pvalue in each permutation
  
  # IMPORTANT: permutations file should have permutations in columns, and column 1 should correspond to the first permutation.
  # specific_permutations should indicate the column/permutation indices of interest such as, 1:10 for permutations 1 to 10
  # is_permutation variable not needed, since the answer is yes.
  # R_obj variable not needed; the R_object for each permutation will always be extracted in order to identify the lowest pvalue.
  # gene variable not needed, since expression data will be defined for each permutation.
  # SNP_file is not needed, since snps will be loaded outside of the main function in order to avoid loading the same snps for each permutation.
  
  
  
  
  
  # LOAD PERMUTATIONS  
  
  # Load file with list of different permutations
  sample_permutations <- read.csv(permutations_file, header = F, row.names = 1, sep = '\t', stringsAsFactors = F)
  
  # Subset specific permutations you want to use and define the total number of permutations
  sample_permutations <- sample_permutations[, specific_permutations, drop = FALSE]
  num_of_perms <- length(specific_permutations)
  
  # Load expression data (only needed to order the permutations to ensure swapping of expression values)
  expression <- read.csv(expression_file, header = T, row.names = 1, sep = '\t', stringsAsFactors = F)
  
  # Re-arrange permutations table so that sample order matches the sample order in the expression file
  sample_permutations <- sample_permutations[colnames(expression), , drop = FALSE]
  
  # Define matrix to hold 1) the lowest pvalue for each permutation 
  lowest_permutations <- matrix(0, nrow = num_of_perms, ncol = 1)
  rownames(lowest_permutations) <- paste("Perm", specific_permutations, sep = "")
  
  
  
  # LOAD GENOTYPE DATA
  # Data is currently loaded outside any function (so it can be used across real/permuted eQTLs)
  
  
  # Run permutation analysis. NOTE: ENSURE 1ST COLUMN IS NOT THE ORIGINAL SAMPLE ORDER
  for (ith_permutation in 1:num_of_perms) {
    
    # Report function progress
    message(paste(ith_permutation-1, 'out of', num_of_perms, 'permutations complete.'))
    
    # Change output file name for each permutation
    ith_R_object_filename <- paste(R_object_filename, '_Perm_', specific_permutations[ith_permutation], sep = "")
    
    # Load expression data (in matrixeqtl format) 
    gene = SlicedData$new()
    gene$fileDelimiter = "\t"      # the TAB character
    gene$fileOmitCharacters = "NA" # denote missing values;
    gene$fileSkipRows = 1          # one row of column labels
    gene$fileSkipColumns = 1       # one column of row labels
    gene$fileSliceSize = 2000      # read file in slices of 2,000 rows
    gene$LoadFile(expression_file)
    
    # Define samples indices in the order that matches the iTH PERMUTATION ORDER 
    expression_indices <- match(sample_permutations[, ith_permutation], gene$columnNames)
    
    # Reorder/subset the samples in the ith PERMUTATION ORDER
    gene$ColumnSubsample(expression_indices)
      
      
      
    # Run matrixeqtleQTL using the randomly swapped expression values
    ith_permutation_eQTL <- Single_Matrix_eQTL(is.perm = TRUE,
                                               output_file_cis = output_file_cis,
                                               snps_pos = snps_pos,
                                               gene_pos = gene_pos,
                                               pvOutputThreshold_cis = pvOutputThreshold_cis,
                                               cis_Dist = cis_Dist,
                                               output_file_tra = output_file_tra,
                                               pvOutputThreshold_tra = pvOutputThreshold_tra,
                                               expression_file = expression_file,
                                               gene = gene,
                                               SNP_file = NULL,
                                               snps = snps,
                                               covariates_file = covariates_file,
                                               R_obj = TRUE, 
                                               output.dir = output.dir,
                                               R_object_filename = ith_R_object_filename)
      
    # Save lowest pvalue for ith permutation
    if ( pvOutputThreshold_tra == 0 ) {
      # this is for cis-eqtl analyses
      
      # Extract lowest pvalue
      lowest_permutations[ith_permutation, 1] <- min(ith_permutation_eQTL$cis$eqtls$pvalue)
        
      # Remove eQTL_object
      rm(ith_permutation_eQTL)

    } else {
      # this is for trans-eqtl analyses
      
      # Extract lowest pvalue
      lowest_permutations[ith_permutation, 1] <- min(ith_permutation_eQTL$all$eqtls$pvalue)
        
      # Remove eQTL_object
      rm(ith_permutation_eQTL)
    }
    

  } # End FOR loop over permutations
    
  
  # Save permutations pvalue vector
  write.table(lowest_permutations, file = (paste(output.dir, R_object_filename, "_lowest_pvalues_", specific_permutations[1], "_", Sys.Date(), ".txt", sep="")), row.names = T, col.names = F, sep = "\t", quote = F)
  
  
  return(NULL)
  
} # END FUNCTION