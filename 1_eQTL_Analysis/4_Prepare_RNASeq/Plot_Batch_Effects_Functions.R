extract_mapping_percents <- function(dir.mapping_files, file.GEUV_metadata){
  # FUNCTION INFO: THIS FUNCTION TAKES ALL THE STAR ALIGNMENT STATISTICS FILES IN THE PROVIDED DIRECTORY AS AN INPUT, AND OUTPUTS A MATRIX WITH UNIQUE/MULTIMAPPING PERCENTS/READ NUMBERS
  # A FILE WITH THE GEUVADIS METADATA FOR EACH SAMPLE IS NEEDED AS AN INPUT
  
  
  # Load sample metadata
  GEUV_sample_info <- read.csv(file.GEUV_metadata, header = T, stringsAsFactors = F, sep = '\t')
  rownames(GEUV_sample_info) <- GEUV_sample_info$Comment.ENA_RUN.
  
  # Set the directory with the STAR logs as the working directory
  setwd(dir.mapping_files)
  
  # List all the STAR log files
  filenames <- Sys.glob("*.final.out")
  
  # Make empty matrix to hold mapping rate values
  mappability <- as.data.frame(matrix(0, nrow = length(filenames), ncol = 6))
  
  # Define counter to keep track of the number of samples
  counter <- 0
  
  
  
  for (file in filenames) {
    
    # Update counter
    counter <- counter + 1
    
    # Read file_i
    file_i_info <- read.csv(file, sep = '|', header = F, row.names = NULL, strip.white = T)
    
    # For row names, use the name of one of the mapping metrics STAR outputs (% of unique reads, % of multi mapping reads, etc.)
    rownames(file_i_info) <- file_i_info$V1
    
    # Fill in the first column with the file name/sample name of the ith file
    mappability[counter, 1] <- filenames[counter]
    
    # Assign various mapping values to columns 2-6, from the ith file
    mappability[counter, c(2,3,4,5,6)] <- file_i_info[c("Uniquely mapped reads %", "% of reads mapped to multiple loci", "Number of input reads", "Uniquely mapped reads number", "Number of reads mapped to multiple loci"), 2]
    
    
  } #End FOR loop
  
  
  
  # Clean up the mappability matrix
  
  # Assign colnames
  colnames(mappability) <- c('Sample', 'Unique_fraction', 'Multimapping_fraction', 'Input_reads', 'Unique_reads', 'Multimapping_reads')
  
  # Remove percentage symbols
  mappability$Unique_fraction <- gsub("%", "", mappability$Unique_fraction, fixed = T)
  mappability$Multimapping_fraction <- gsub("%", "", mappability$Multimapping_fraction, fixed = T)
  
  # Update the percent columns to numeric values, and convert percents to fractions
  mappability$Unique_fraction <- as.numeric(mappability$Unique_fraction)/100
  mappability$Multimapping_fraction <- as.numeric(mappability$Multimapping_fraction)/100
  
  # Shorten sample names
  mappability$Sample <- sub("Log.final.out*", "", mappability$Sample)
  mappability$Sample <- sub("fastp_*", "", mappability$Sample)
  
  # Get alternative samples names used by Geuvadis
  rownames(mappability) <- GEUV_sample_info[mappability$Sample, 'Source.Name']
  


  
  
  # Return output
  return(mappability)
  
  
} # END FUNCTION  