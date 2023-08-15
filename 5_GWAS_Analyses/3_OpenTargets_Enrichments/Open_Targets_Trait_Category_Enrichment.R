# Set strings as factors
options(stringsAsFactors = F)

# Load required libraries
library(arrow)
library(data.table)
library(tidyverse)
library(dplyr)
library(ghql)
library(jsonlite)
library(magrittr)
library(rlang)

# Load functions
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/5_GWAS_Analyses/3_OpenTargets_Enrichments/Open_Targets_Trait_Category_Enrichment_functions.R")

# Specify the directory for parquet mapping files
dir.parquet <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/5_GWAS_Analyses/3_OpenTargets_Enrichments/RsID_to_OpenTargets_Maps/'

# Specify directory for results
dir.results <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/5_GWAS_Analyses/3_OpenTargets_Enrichments/Results/'





    
    
    
    
  
# GENERATE A DATAFRAME OF RSIDS IN OPENTARGETS
    
# Download parquet files via termina
# wget --recursive --no-parent --no-host-directories --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/lut/variant-index/ 

# Collect individual file paths in a list.
parquet_list <- list.files(dir.parquet, "\\.parquet$", recursive=FALSE, full.names=TRUE)
  
# Create new dataframe to hold all Open Targets SNPs
OpenTargets_variants <- data.frame()

# Start a counter
parquet_counter <- 0

# Loop over each parquet file path
for(parquet_file in parquet_list) {
  
    # Update the parquet counter
    parquet_counter <- parquet_counter + 1
    
    # Read in the ith parquet file
    ith_parquet_data <- read_parquet(parquet_file)
    
    # Filter SNPs that are not bi-allelic
    
        # Count the number of bases for each REF allele
        ith_parquet_data$REF_length <- apply(ith_parquet_data[, "ref_allele"], 2, nchar)
        
        # Count the number of bases for each ALT allele
        ith_parquet_data$ALT_length <- apply(ith_parquet_data[, "alt_allele"], 2, nchar)
        
        # Only keep OpenTarget SNPs with 1 base for each of the REF/ALT alleles (i.e. remove entries with indels)
        ith_parquet_data <- ith_parquet_data[which(ith_parquet_data$REF_length == 1), ]
        ith_parquet_data <- ith_parquet_data[which(ith_parquet_data$ALT_length == 1), ]
    
    # Define or extend a df with ith data
    if (parquet_counter == 1) {
      
      # If this is the first parquet file, define a new df to hold the data. Keep only the first 12 columns.
      OpenTargets_variants <- ith_parquet_data[, 1:12]
      
    } else {
      
      # If this is NOT the first parquet file, add the data to the existing df. Keep only the first 12 columns.
      OpenTargets_variants <- rbind(OpenTargets_variants, ith_parquet_data[, 1:12])
      
    }
  
} # END FOR LOOP

# Remove temporary variables
rm(ith_parquet_data)

# Save the csv
fwrite(OpenTargets_variants, paste(dir.parquet, "Open_Targets_Variant_Info.csv", sep = ''))










# DEFINE OPEN TARGETS IDS FOR ALL RSIDS (BOTH BACKGROUND AND SIGNIFICANT SNPS)

# Read in all SNVs tested in the eQTL analysis
all_snps <- fread(file="/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/RESOURCE_SNP_info_EUR.txt", header = TRUE, sep = '\t')

    # Update colnames
    colnames(all_snps)[which(names(all_snps) == "RsID")] <- "rs_id"
    
# Read in significant trans-eQTL SNPs found by the lab
sig_snps <- read.csv(file="/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/5_GWAS_Analyses/1_Prepare_Files/Bravo_eQTL_SNVs_for_Collaborators_2023-05-30.txt", header=TRUE, sep='\t') 

    # Update colnames
    colnames(sig_snps)[which(names(sig_snps) == "SNP")] <- "RsID"
    colnames(sig_snps)[which(names(sig_snps) == "HG38_CHR")] <- "chr"
    colnames(sig_snps)[which(names(sig_snps) == "HG38_BP")] <- "pos"
    colnames(sig_snps)[which(names(sig_snps) == "RsID")] <- "rsid"
    
# Extract Open Targets data for all background eQTL SNPs    
all_SNVs_in_open_targets <- semi_join(OpenTargets_variants, all_snps, by = "rs_id")

# Update the RsID colname
colnames(all_SNVs_in_open_targets)[which(names(all_SNVs_in_open_targets) == "rs_id")] <- "rsid"

# Define a column with SNP ID in Open Targets format
all_SNVs_in_open_targets$variant_id <- paste(all_SNVs_in_open_targets$chr_id, 
                                             all_SNVs_in_open_targets$position, 
                                             all_SNVs_in_open_targets$ref_allele, 
                                             all_SNVs_in_open_targets$alt_allele, sep="_") 

# Keep only the ID columns
all_SNVs_in_open_targets <- all_SNVs_in_open_targets[, c('rsid', 'variant_id')]

    # Save as a csv file
    fwrite(all_SNVs_in_open_targets, paste(dir.parquet, "Open_Targets_IDs_All_eQTL_SNVs.csv", sep = ''))


# Extract Open Targets IDs for all significant eQTL SNPs 
significant_SNVs_in_open_targets <- semi_join(all_SNVs_in_open_targets, sig_snps, by = "rsid")

    # Save as a csv file
    fwrite(significant_SNVs_in_open_targets, paste(dir.parquet, "Open_Targets_IDs_Significant_eQTL_SNVs.csv", sep = ''))

# Remove temporary variables
rm(OpenTargets_variants, all_snps, sig_snps)










# GENERATE RANDOM SNV SAMPLES TO EMPIRICALLY CALCULATE TRAIT CATEGORY ENRICHMENTS

# Read in background SNV OpenTargets IDs
all_SNVs_in_open_targets <- fread(file="/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/5_GWAS_Analyses/3_OpenTargets_Enrichments/RsID_to_OpenTargets_Maps/Open_Targets_IDs_All_eQTL_SNVs.csv", header = TRUE, sep = ',')

# Read in significant SNV OpenTargets IDs
significant_SNVs_in_open_targets <- fread(file="/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/5_GWAS_Analyses/3_OpenTargets_Enrichments/RsID_to_OpenTargets_Maps/Open_Targets_IDs_Significant_eQTL_SNVs.csv", header = TRUE, sep = ',')

# Set the number of random combinations we want to generate
num_random_samples <- 500 

# Set the number of SNVs in each sample
numb_per_sample <- nrow(significant_SNVs_in_open_targets)

# Create a dataframe to hold combinations. Rows are combinations. Generate 2X combinations than necessary, in case there are duplicates.
my.random.combinations <- data.frame(matrix(nrow = 2 * num_random_samples, ncol = numb_per_sample))

# Set random seed
set.seed(12345) 

# Generate combinations 
for(ith_combination in 1:nrow(my.random.combinations)){
  
  # Randomly sample indices for the ith combination
  row_indices <- sample.int(n = nrow(all_SNVs_in_open_targets), size = numb_per_sample, replace = FALSE)
  
  # Order indices from smallest to largest
  row_indices <- sort(row_indices, decreasing = FALSE)
  
  # Add the combination to the df
  my.random.combinations[ith_combination, ] <- t(all_SNVs_in_open_targets[row_indices, "variant_id"]) 
  
  # Print the ith combination
  print(ith_combination)
}

# Remove duplicate rows
my.random.combinations <- my.random.combinations[!duplicated(my.random.combinations), ]

# Keep only the first N unique combinations
my.random.combinations <- my.random.combinations[1:num_random_samples, ]

# Save as a csv file
fwrite(my.random.combinations, paste(dir.parquet, "Open_Targets_ID_Random_Combinations.csv", sep = ''))










# OBTAIN PHENOTYPE ASSOCIATIONS FOR L1 EQTLS

# Read in significant SNV OpenTargets IDs
significant_SNVs_in_open_targets <- fread(file="/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/5_GWAS_Analyses/3_OpenTargets_Enrichments/RsID_to_OpenTargets_Maps/Open_Targets_IDs_Significant_eQTL_SNVs.csv", header = TRUE, sep = ',')

# Call the open targets graphql API
otg_cli <- GraphqlClient$new(url = "https://api.genetics.opentargets.org/graphql") 

    # Create a new query
    otg_qry <- Query$new() 

# Get significant pheWAS studies and traits from the significant L1 eQTLs
my_sig_phewas <- get_phewas(significant_SNVs_in_open_targets, pval_filter = 5e-8) 

    # Save as a csv file
    fwrite(my_sig_phewas, paste(dir.results, "All_Significant_PheWAS_Results_for_L1_eQTLs.csv", sep = ''))

# Define trait categories we want to test for enrichment
list_trait_categories <- c("endocrine system disease", "musculoskeletal or connective tissue disease", "integumentary system disease", "gastrointestinal disease", "urinary system disease",
                           "disease of visual system", "hematologic disease", "immune system disease", "pancreas disease", "cell proliferation disorder", "nervous system disease",
                           "respiratory or thoracic disease", "cardiovascular disease", "disease of ear")

# Define trait categories that won't be included (for reference purposes only)
list_omitted_trait_categories <- c("Uncategorised", "pregnancy or perinatal disease", "phenotype", "measurement", "injury, poisoning or other complication", "infectious disease", "congenital disease", "biological process")

# In the phewas results, only keep unique combinations of rsid and traitCategory
my_sig_phewas_categories <- my_sig_phewas %>% distinct(rsid, traitCategory, .keep_all=TRUE)

# Create a dataframe that will include the number of SNPs mapping to each trait category
phewas_category_counts_eqtls <- data.frame(matrix(ncol = length(list_trait_categories), nrow = 1)) 
colnames(phewas_category_counts_eqtls) <- setNames(as.list(list_trait_categories), list_trait_categories) 

# Count the number of SNPs in the eQTL list that are associated with each category
for (ith_category in list_trait_categories) {
  
    # Count the number of SNPs mapping to the ith category
    ith_category_count <- sum(my_sig_phewas_categories$traitCategory == ith_category)
    
    # Add the count to the count table
    phewas_category_counts_eqtls[1, ith_category] <- ith_category_count
  
}

# Save as a csv file
fwrite(phewas_category_counts_eqtls, paste(dir.results, "Open_Targets_Disease_Categories_eQTL_SNP_Counts.csv", sep = ''))










# OBTAIN PHENOTYPE ASSOCIATIONS FOR RANDOM COMBINATIONS

# Read in disease category counts for eQTLs
phewas_category_counts_eqtls <- fread(file="/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/5_GWAS_Analyses/3_OpenTargets_Enrichments/Results/Open_Targets_Disease_Categories_eQTL_SNP_Counts.csv", header = TRUE, sep = ',')

# Read in file with random combinations of SNPs
my.random.combinations <- fread(file="/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/5_GWAS_Analyses/3_OpenTargets_Enrichments/RsID_to_OpenTargets_Maps/Open_Targets_ID_Random_Combinations.csv", header = TRUE, sep = ',')

# Create an empty data frame to hold category counts for the random combinations
phewas_category_counts_RANDOM <- data.frame(matrix(ncol = ncol(phewas_category_counts_eqtls), nrow = nrow(my.random.combinations)))

    # Add the column names of the trait categories, the rows will correspond with the counts of that category for each random combination
    colnames(phewas_category_counts_RANDOM) <- colnames(phewas_category_counts_eqtls)
    rownames(phewas_category_counts_RANDOM) <- paste("Random_Combination_", 1:nrow(phewas_category_counts_RANDOM), sep = '')

# Call the open targets graphql API
otg_cli <- GraphqlClient$new(url = "https://api.genetics.opentargets.org/graphql") 

    # Create a new query
    otg_qry <- Query$new() 
    
# Loop over each random combination
for(ith_combination in 1:nrow(my.random.combinations)){ 
  
    # Define the start time
    start_time <- Sys.time()
    
    # Get significant PheWAS results for the ith combination
    ith_phewas <- get_phewas(my.random.combinations[ith_combination, ], pval_filter = 5e-8)
    
    # Only keep unique combinations of variantInfo and traitCategory
    ith_phewas_categories <- ith_phewas %>% distinct(variantInfo, traitCategory, .keep_all = TRUE)

    # For each category in list_trait_categories
    for (ith_category in colnames(phewas_category_counts_RANDOM)) {

        # Count the number of times each category appears in the ith combination
        ith_random_category_count <- sum(ith_phewas_categories$traitCategory == ith_category)

        # Add the count to the right trait category
        phewas_category_counts_RANDOM[ith_combination, ith_category] <- ith_random_category_count

    }

    # Define the end time
    end_time <- Sys.time()

    # Define the elapsed time
    elapsed_time <- end_time - start_time

    # Print the elapsed time
    print(paste("Elapsed time for iteration", ith_combination, ":", elapsed_time))


} # END FOR LOOP OVER COMBINATIONS

# Save as a csv file
fwrite(phewas_category_counts_RANDOM, paste(dir.results, "Open_Targets_Disease_Categories_Random_SNP_Counts.csv", sep = ''), row.names = TRUE)







# EMPIRICALLY CALCULATE DISEASE CATEGORY ENRICHMENTS AND GENERATE PLOTS

# Read in disease category counts for eQTLs
phewas_category_counts_eqtls <- data.frame(fread(file="/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/5_GWAS_Analyses/3_OpenTargets_Enrichments/Results/Open_Targets_Disease_Categories_eQTL_SNP_Counts.csv", header = TRUE, sep = ','))

# Read in disease category counts for combinations
phewas_category_counts_RANDOM <- data.frame(fread(file="/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/5_GWAS_Analyses/3_OpenTargets_Enrichments/Results/Open_Targets_Disease_Categories_Random_SNP_Counts.csv", header = TRUE, sep = ','), row.names = 1)

# Define a df to hold the stats
Enrichment_Stats <- data.frame(matrix(nrow = ncol(phewas_category_counts_eqtls), ncol = 2))

# Update rownames and colnames
colnames(Enrichment_Stats) <- c("p_value", "FDR")
rownames(Enrichment_Stats) <- colnames(phewas_category_counts_RANDOM)

# For each category, find the p-value for the significant SNPs
for(ith_category in colnames(phewas_category_counts_RANDOM)) {
  
  # Extract random data for the ith category
  ith_category_random_data <- phewas_category_counts_RANDOM[, ith_category]
  
  # Calculate the pvalue for the ith category
  ith_p_value <- get_p_value(real_data_count = phewas_category_counts_eqtls[1, ith_category], 
                             random_data_counts = ith_category_random_data,
                             num_perms = nrow(phewas_category_counts_RANDOM)) 
  
  # Add the pvalue to the enrichment stats table
  Enrichment_Stats[ith_category, "p_value"] <- ith_p_value
}

# Adjust p-values
Enrichment_Stats[, "FDR"] <- p.adjust(Enrichment_Stats$p_value, method = 'fdr')

    # Save as a csv file
    fwrite(Enrichment_Stats, paste(dir.results, "Open_Targets_Disease_Categories_Enrichment_Stats.csv", sep = ''), row.names = TRUE)
    
    

# Merge the real and random data (to use in selecting plot bounds)
combined_real_random <- rbind(phewas_category_counts_RANDOM, phewas_category_counts_eqtls)

# For each category, create an ecdf plot using the samples and plot where the significant SNP values fall along with FDR value
for (ith_category in colnames(phewas_category_counts_RANDOM)) {
  
  # Define file name
  file <- paste(dir.results, "ECDF_Plot_", ith_category, ".pdf", sep = "")
  
  # Generate a pdf
  pdf(file, width = 5, height = 5)
  
   # Extract the eqtl count for the ith category
  ith_category_real_value <- phewas_category_counts_eqtls[1, ith_category]
  
  # Extract random data for the ith category
  ith_category_random_data <- phewas_category_counts_RANDOM[, ith_category]
  
  # Order random data counts from smallest to largest
  ith_category_random_data <- sort(ith_category_random_data, decreasing = FALSE)
  
  # Generate a df to hold manually calculated cumulative probabilities
  manual_cumulative <- data_frame(counts = ith_category_random_data, cumulative_p = 0)
  
  # Calculate and fill in cumulative pvalues
  for (ith_combination in 1:length(ith_category_random_data)) {
    
    # Update the cumulative pvalue
    manual_cumulative[ith_combination, "cumulative_p"] <- ith_combination/length(ith_category_random_data)
    
  }
  
  # Generate an ecdf function for the ith category
  my_ith_ecdf <- ecdf(ith_category_random_data)
  
  # Calculate cumulative probability for the ith category real data count
  ith_category_cumulative_p <- my_ith_ecdf(ith_category_real_value)
  
  # Specify the FDR, limiting the value to 2 significant figures
  ith_FDR <- Enrichment_Stats[ith_category, "FDR"]
  ith_FDR <- formatC(ith_FDR, format = 'E', digits = 2)
  
  # Calculate an enrichment score (real data counts/median of the random combinations)
  ith_enrichment <- ith_category_real_value/median(ith_category_random_data)
  
  # Generate a plot
  plot(manual_cumulative, 
       main = ith_category, 
       xlab = "Number of SNPs Mapping to Disease Category", 
       ylab = "Cumulative Probability", 
       xlim = c(0, max(combined_real_random[, ith_category])+30), 
       ylim =c (0, 1.1),
       pch = 20, 
       cex = 0.5)
  
  # overlay the point for the real data point
  points(x = ith_category_real_value, y = ith_category_cumulative_p, col = 'red', pch = 20, cex = 0.5)

  # Draw a horizontal line at y = 0 and at y = 1
  abline(h = 0, col = "grey", lty = "dashed", lwd = 1)
  abline(h = 1, col = "grey", lty = "dashed", lwd = 1)
  
  # Add text specifying real data count, FDR, and enrichment score
  text(x = ith_category_real_value, y = ith_category_cumulative_p-0.1, labels = paste("Value:", ith_category_real_value, "\nFDR:", ith_FDR, "\nES:", ith_enrichment), pos = 4, col = "red")
  
  # End plot
  dev.off()
  
}



# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/5_GWAS_Analyses/3_OpenTargets_Enrichments/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_OpenTargets_Enrichment.txt", sep =""))
sessionInfo()
sink()      
    



