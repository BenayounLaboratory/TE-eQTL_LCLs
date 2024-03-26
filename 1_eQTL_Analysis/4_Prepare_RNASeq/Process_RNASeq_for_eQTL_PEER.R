# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Process_RNASeq_for_eQTL_functions.R")

# Load libraries
#library(edgeR)
#library(genefilter)

#library(DESeq2) # For VST
#library(BiocParallel) # Used by DESeq2 for parallelization
#  register(MulticoreParam(6)) # Use six cores
library(limma) # for batch effect removal
library(RNOmni) # for inverse normal transform

  
# Define output directory
dir.output <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_with_PEER/'

    
    
    
# Section 1: INFER HIDDEN COVARIATES WITH PEER (OPEN DOCKER AND RUN THIS IN TERMINAL) ------------------------------------------------------------------------------------------------------------------------
 

# Setup PEER Docker

  # Download the docker image
  # Docker pull quay.io/biocontainers/r-peer:1.3--r341h470a237_1  

  # Mount the volume into my_data container and run the docker image
  # Docker run -it --name test2 -v /Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_with_PEER:/Processed_counts_with_PEER quay.io/biocontainers/r-peer:1.3--r341h470a237_1

  # Run R in docker image (R is verion 3.4.1)
  # /usr/local/lib/R/bin/R

  
## NOTE: THIS NEEDS TO BE RUN ON DOCKER IMAGE WITH PEER, SINCE IT IS NOT COMPATIBLE WITH CURRENT R VERSION
  
# Load PEER library
library(peer)

# Set the working directory
setwd('/Processed_counts_with_PEER/')

# Load the expression data
expr.EUR = read.csv('Inputs/All_counts_EUR_358_filtered_VST.txt', header = TRUE, row.names = 1, sep = '\t')

    # Transpose the data (rows should be samples and columns should be genes)
    expr.EUR <- t(expr.EUR)

    # Check expression data dimensions
    dim(expr.EUR) 
    
# Initiate a PEER model
my.model = PEER()

# Assign the expression data to the model
PEER_setPhenoMean(my.model, as.matrix(expr.EUR)) # NULL means no errors

# Read-in known covariates
covs.EUR = read.csv('Inputs/COVARIATES_EUR.txt', header = TRUE, row.names = 1, sep = '\t')

    # Ensure samples are in the same order as expression
    covs.EUR <- covs.EUR[rownames(expr.EUR), ]  

    # Check dimensions and column names
    dim(covs.EUR) 
    colnames(covs.EUR)
    
    # Remove the first column (the sample ID column)
    covs.EUR <- covs.EUR[, -c(1)]
    
    # Change categorical variables to binary variables
    binary.cov.EUR <- model.matrix(~ sex + ancestry + lab + PC1 + PC2 + Net_L1_Alu_Copies + EBV_expr_VST, data = covs.EUR)
    
    # Remove the first column (intercept)
    binary.cov.EUR <- binary.cov.EUR[, -c(1)]

# Assign the covariates to the model
PEER_setCovariates(my.model, as.matrix(binary.cov.EUR))

# Define the number of hidden confounders to infer.The GEUVADIS paper didn't use known covariates, only 10 PEER factors.
PEER_setNk(my.model, 10) 

# Perform the inference (219 iterations for convergence)
PEER_update(my.model)

# Obtain outputs
my.factors = PEER_getX(my.model)
dim(my.factors)
write.table(my.factors, file = (paste("PEER_output/PEER_FACTORS", ".txt", sep="")), row.names = F, col.names = F, sep = "\t", quote = F)

my.residuals = PEER_getResiduals(my.model)
dim(my.residuals)
write.table(my.residuals, file = (paste("PEER_output/PEER_RESIDUALS", ".txt", sep="")), row.names = F, col.names = F, sep = "\t", quote = F)

my.weights = PEER_getW(my.model)
dim(my.weights)
write.table(my.weights, file = (paste("PEER_output/PEER_WEIGHTS", ".txt", sep="")), row.names = F, col.names = F, sep = "\t", quote = F)

my.precision = PEER_getAlpha(my.model)
dim(my.precision)
write.table(my.precision, file = (paste("PEER_output/PEER_PRECISION", ".txt", sep="")), row.names = F, col.names = F, sep = "\t", quote = F)

# Generate diagnostic plots
pdf('PEER_output/PEER_Plots_1.pdf')
PEER_plotModel(my.model)
dev.off()

pdf('PEER_output/PEER_Plots_2.pdf')
plot(my.precision)
dev.off()
    


  

# Section 2: BATCH EFFECT REMOVAL WITH LIMMA -------------------------------------


# BATCH EFFECT REMOVAL WITH LIMMA

# Load VST expression data
VST.EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_eQTL/All_counts_EUR_358_filtered_VST.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Load the covariates data
covariates_table.EUR  <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_eQTL/COVARIATES_EUR.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
covariates_table.EUR <- covariates_table.EUR[colnames(VST.EUR), ]   
    
    # Make model matrices
    full.model.EUR = model.matrix(~ sex + ancestry + lab + PC1 + PC2 + Net_L1_Alu_Copies + EBV_expr_VST, data = covariates_table.EUR)

# Load PEER factors, subset PEER columns, and assign column names.
PEER_EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_with_PEER/PEER_output/PEER_FACTORS.txt", header = F, stringsAsFactors = F, sep = '\t')
PEER_EUR <- PEER_EUR[tail(names(PEER_EUR), 10)] # 10 PEER Covariates. IF THIS IS CHANGED,MAKE SURE IT IS CORRECT
colnames(PEER_EUR) <- paste('PEER', 1:length(PEER_EUR), sep = '')

# Add PEER factors to covariates table
full.model.EUR <- cbind(full.model.EUR, PEER_EUR) 

# From the VST data, regress out lab + ancestry + PCs + sex + TE_Insertions + EBV + PEER
batch_corr_EUR <- removeBatchEffect(x = VST.EUR,
                                    batch = NULL,
                                    covariates = full.model.EUR[, c(2:25)],
                                    design = full.model.EUR[, c(1), drop = FALSE])

# Save batch corrected data
write.table(batch_corr_EUR, file = paste(dir.output, 'All_counts_EUR_358_filtered_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV_PEER', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)


            

            
# Section 3: INVERSE NORMAL TRANSFORMATIONS (INT) -------------------------------------
# Inverse Normal Transformation (INT), so data is normally distributed, as needed for linear regression in eqtl
        

# Load batch-corrected VST expression data
batch_corr_EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_with_PEER/All_counts_EUR_358_filtered_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV_PEER.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Run function to calculate INT expression for each gene (THESE VALUES WILL SERVE AS THE INPUT FOR MATRIXEQTL)
INT.vst.regressed.EUR <- Inverse_normal_transform(input.df = batch_corr_EUR)

# Run function to save INT matrices, and also generate subsets by gene, L1 subfamilies, and Alu subfamilies
Split_and_save_expression(expression_df = INT.vst.regressed.EUR, 
                          output.dir = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_with_PEER/', 
                          output_file_prefix = 'All_counts_EUR_358_filtered_VST_BatchesRemoved_lab_ancestry_PCs_sex_TE_copies_EBV_PEER_INT_', 
                          gene_prefix = 'ENSG', 
                          TE_groups = 'Subfamily')


          
# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Session_Info/'
    
    
sink(file = paste(dir.session_info,"Session_Info_Process_RNASeq_for_PEER.txt", sep =""))
sessionInfo()
sink()    
    


# Clean the environment
rm(list=ls())


