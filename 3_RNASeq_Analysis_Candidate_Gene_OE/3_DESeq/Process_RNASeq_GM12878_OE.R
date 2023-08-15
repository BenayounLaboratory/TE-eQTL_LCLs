# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/3_DESeq/Process_RNASeq_functions.R")

# Load libraries
library(DESeq2) # For differential expression
library(BiocParallel) # Used by DESeq2 for parallelization
  register(MulticoreParam(6)) # Use six cores
library(sva) # Correcting for batch effects
library(limma) # Correcting for batch effects
library(pheatmap) # For gene expression heatmaps
library(biomaRt) # To map Ensembl IDs to gene symbols


    
    
    
# Section 1: FILTER LOW EXPRESSION GENES ------------------------------------------------------------------------------------------------------------------------
 


  
# Define the output directory
dir.output <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/3_DESeq/Processed_counts_GM12878_OE/'

# Define the directory with the counts tables and collect individual file paths in a list.
count_Table.dir <- c("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/2_Read_counting/Counts_GM12878_Candidate_OE/")
cntTable_list <- list.files(count_Table.dir, "\\.cntTable$", recursive=FALSE, full.names=TRUE)


# Run function to combine the individual counts tables into one df.
gene_TE_counts_raw <- aggregate_counts(cntTable_complete.list = cntTable_list)


# Remove useless info. from sample names
colnames(gene_TE_counts_raw) <- sub("Aligned.sortedByCoord.out.bam*", "", colnames(gene_TE_counts_raw))
colnames(gene_TE_counts_raw) <- sub("*fastp_", "", colnames(gene_TE_counts_raw))


# Reorder columns so controls are first
gene_TE_counts_raw <- gene_TE_counts_raw[, c('gene.TE', 
                                             'CAN1', 'CAN6', 'CAN11', 'CAN16',
                                             'CAN2', 'CAN7', 'CAN12', 'CAN17',
                                             'CAN3', 'CAN8', 'CAN13', 'CAN18',
                                             'CAN4', 'CAN9', 'CAN14', 'CAN19',
                                             'CAN5', 'CAN10', 'CAN15', 'CAN20')]
      
## Assign more specific sample names
colnames(gene_TE_counts_raw) <- c('gene.TE', 
                                  'Empty_1', 'Empty_2', 'Empty_3', 'Empty_4',
                                  'IL16_1', 'IL16_2', 'IL16_3', 'IL16_4',
                                  'STARD5_1', 'STARD5_2', 'STARD5_3', 'STARD5_4',
                                  'HSD17B12_1', 'HSD17B12_2', 'HSD17B12_3', 'HSD17B12_4',
                                  'RNF5_1', 'RNF5_2', 'RNF5_3', 'RNF5_4')

# Remove non-gene/TE features from expression table
misc.feature.indices <- which(!grepl('EBV', gene_TE_counts_raw$gene.TE) & !grepl('pcDNA', gene_TE_counts_raw$gene.TE))
gene_TE_counts_raw <- gene_TE_counts_raw[misc.feature.indices, ]


# Define the controls/treatments for each analysis
treatment1 <- c('gene.TE', 
                'Empty_1', 'Empty_2', 'Empty_3', 'Empty_4',
                'IL16_1', 'IL16_2', 'IL16_3', 'IL16_4')

treatment2 <- c('gene.TE', 
                'Empty_1', 'Empty_2', 'Empty_3', 'Empty_4',
                'STARD5_1', 'STARD5_2', 'STARD5_3', 'STARD5_4'
               )

treatment3 <- c('gene.TE', 
                'Empty_1', 'Empty_2', 'Empty_3', 'Empty_4',
                'HSD17B12_1', 'HSD17B12_2', 'HSD17B12_3', 'HSD17B12_4')

treatment4 <- c('gene.TE', 
                'Empty_1', 'Empty_2', 'Empty_3', 'Empty_4',
                'RNF5_1', 'RNF5_2', 'RNF5_3', 'RNF5_4')

# Run function to remove gene/transcript version info, combine transcripts from the same gene, and filter lowly expressed genes
counts_filtered.IL16 <- cleanup_and_filter_counts(count_data = gene_TE_counts_raw[, treatment1], filter = TRUE, min_counts_per_sample = 10, fraction_of_samples = 4/8)
counts_filtered.STARD5 <- cleanup_and_filter_counts(count_data = gene_TE_counts_raw[, treatment2], filter = TRUE, min_counts_per_sample = 10, fraction_of_samples = 4/8)
counts_filtered.HSD17B12 <- cleanup_and_filter_counts(count_data = gene_TE_counts_raw[, treatment3], filter = TRUE, min_counts_per_sample = 10, fraction_of_samples = 4/8)
counts_filtered.RNF5 <- cleanup_and_filter_counts(count_data = gene_TE_counts_raw[, treatment4], filter = TRUE, min_counts_per_sample = 10, fraction_of_samples = 4/8)


    # Save counts files.
    write.table(counts_filtered.IL16, file = paste(dir.output, "All_counts_filtered_IL16", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)
    write.table(counts_filtered.STARD5, file = paste(dir.output, "All_counts_filtered_STARD5", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)
    write.table(counts_filtered.HSD17B12, file = paste(dir.output, "All_counts_filtered_HSD17B12", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)
    write.table(counts_filtered.RNF5, file = paste(dir.output, "All_counts_filtered_RNF5", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)


    
    
    


    






# Section 2: SVA + DESeq + MDS/PCA ------------------------------------------------------------------------------------------------------------------------




# Define the output directories
counts.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/3_DESeq/Processed_counts_GM12878_OE/'
DESeq.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/3_DESeq/DESeq_Results_GM12878_OE/'
MDS.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/3_DESeq/MDS_PCA_GM12878_OE/'



# Run SVA for each pairwise analysis

    # IL16
    counts_SVA.IL16 <- Run_SVA(filtered_counts = counts_filtered.IL16, 
                                    control_label = 'Empty', 
                                    control_reps = 4, 
                                    treatment_label = 'IL16', 
                                    treatment_reps = 4, 
                                    output_location = counts.output.dir)
    
    # STARD5
    counts_SVA.STARD5 <- Run_SVA(filtered_counts = counts_filtered.STARD5, 
                                      control_label = 'Empty', 
                                      control_reps = 4, 
                                      treatment_label = 'STARD5', 
                                      treatment_reps = 4, 
                                      output_location = counts.output.dir)
    
    # HSD17B12
    counts_SVA.HSD17B12 <- Run_SVA(filtered_counts = counts_filtered.HSD17B12, 
                                        control_label = 'Empty', 
                                        control_reps = 4, 
                                        treatment_label = 'HSD17B12', 
                                        treatment_reps = 4, 
                                        output_location = counts.output.dir)
    
    # RNF5
    counts_SVA.RNF5 <- Run_SVA(filtered_counts = counts_filtered.RNF5, 
                                    control_label = 'Empty', 
                                    control_reps = 4, 
                                    treatment_label = 'RNF5', 
                                    treatment_reps = 4, 
                                    output_location = counts.output.dir)
    
    
    
# Run DESeq for each treatment; output the VST data and significant genes in a list

    # IL16
    DESeq.IL16 <- Run_DESeq(filtered_counts = counts_SVA.IL16, 
                          padj_limit = 0.05, 
                          control_label = 'Empty', 
                          control_reps = 4, 
                          treatment_label = 'IL16', 
                          treatment_reps = 4, 
                          VST_output = counts.output.dir, 
                          DESeq_output = DESeq.output.dir)
    
    # STARD5
    DESeq.STARD5 <- Run_DESeq(filtered_counts = counts_SVA.STARD5, 
                            padj_limit = 0.05, 
                            control_label = 'Empty', 
                            control_reps = 4, 
                            treatment_label = 'STARD5', 
                            treatment_reps = 4, 
                            VST_output = counts.output.dir, 
                            DESeq_output = DESeq.output.dir)
    
    # HSD17B12
    DESeq.HSD17B12 <- Run_DESeq(filtered_counts = counts_SVA.HSD17B12, 
                              padj_limit = 0.05, 
                              control_label = 'Empty', 
                              control_reps = 4, 
                              treatment_label = 'HSD17B12', 
                              treatment_reps = 4, 
                              VST_output = counts.output.dir, 
                              DESeq_output = DESeq.output.dir)
    
    # RNF5
    DESeq.RNF5 <- Run_DESeq(filtered_counts = counts_SVA.RNF5, 
                          padj_limit = 0.05, 
                          control_label = 'Empty', 
                          control_reps = 4, 
                          treatment_label = 'RNF5', 
                          treatment_reps = 4, 
                          VST_output = counts.output.dir, 
                          DESeq_output = DESeq.output.dir)


# Generate expression box plots for my targets genes
    
    
    # IL16
    Plot_Expression(output.dir = DESeq.output.dir, 
                    target_gene_Ensembl_name = 'ENSG00000172349', 
                    target_gene_symbol = 'IL16', 
                    VST.data = DESeq.IL16[[1]], 
                    Sig.DESeq.Res = DESeq.IL16[[2]],
                    control_label = 'Empty', 
                    control_n = '4', 
                    treatment_label = 'IL16 OE', 
                    treatment_n = '4',
                    ymin = 10.5,
                    ymax = 12.5)
    
    # STARD5
    Plot_Expression(output.dir = DESeq.output.dir, 
                  target_gene_Ensembl_name = 'ENSG00000172345', 
                  target_gene_symbol = 'STARD5', 
                  VST.data = DESeq.STARD5[[1]], 
                  Sig.DESeq.Res = DESeq.STARD5[[2]],
                  control_label = 'Empty', 
                  control_n = '4', 
                  treatment_label = 'STARD5 OE', 
                  treatment_n = '4',
                  ymin = 8.5,
                  ymax = 11)
    
    
    # HSD17B12
    Plot_Expression(output.dir = DESeq.output.dir, 
                  target_gene_Ensembl_name = 'ENSG00000149084', 
                  target_gene_symbol = 'HSD17B12', 
                  VST.data = DESeq.HSD17B12[[1]], 
                  Sig.DESeq.Res = DESeq.HSD17B12[[2]],
                  control_label = 'Empty', 
                  control_n = '4', 
                  treatment_label = 'HSD17B12 OE', 
                  treatment_n = '4',
                  ymin = 10,
                  ymax = 11)
    
    # RNF5
    Plot_Expression(output.dir = DESeq.output.dir, 
                target_gene_Ensembl_name = 'ENSG00000204308', 
                target_gene_symbol = 'RNF5', 
                VST.data = DESeq.RNF5[[1]], 
                Sig.DESeq.Res = DESeq.RNF5[[2]],
                control_label = 'Empty', 
                control_n = '4', 
                treatment_label = 'RNF5 OE', 
                treatment_n = '4',
                ymin = 10,
                ymax = 11.5)
    
    
# Generate express heatmaps for all candidate genes

    # Define genes to plot
    my.candidates <- c('ENSG00000172349',
                       'ENSG00000172345',
                       'ENSG00000149084',
                       'ENSG00000204308')
    
    # Define genes to plot as symbols
    my.candidates.symbols <- c('IL16',
                               'STARD5',
                               'HSD17B12',
                               'RNF5')
    
    # IL16 OE
    Candidate_gene_heatmaps(output.dir = DESeq.output.dir, 
                            EnsemblIDs_to_plot = my.candidates, 
                            VST.data = DESeq.IL16[[1]],
                            my.gene.symbols = my.candidates.symbols, 
                            control_label = 'Empty', 
                            control_n = 4, 
                            control_color = 'black', 
                            treatment_label = 'IL16_OE', 
                            treatment_n = 4, 
                            treatment_color = 'red1')
    
    # STARD5 OE
    Candidate_gene_heatmaps(output.dir = DESeq.output.dir, 
                            EnsemblIDs_to_plot = my.candidates, 
                            VST.data = DESeq.STARD5[[1]],
                            my.gene.symbols = my.candidates.symbols, 
                            control_label = 'Empty', 
                            control_n = 4, 
                            control_color = 'black', 
                            treatment_label = 'STARD5_OE', 
                            treatment_n = 4, 
                            treatment_color = 'pink')
    
    # HSD17B12 OE
    Candidate_gene_heatmaps(output.dir = DESeq.output.dir, 
                            EnsemblIDs_to_plot = my.candidates, 
                            VST.data = DESeq.HSD17B12[[1]],
                            my.gene.symbols = my.candidates.symbols, 
                            control_label = 'Empty', 
                            control_n = 4, 
                            control_color = 'black', 
                            treatment_label = 'HSD17B12_OE', 
                            treatment_n = 4, 
                            treatment_color = '#369193')
    
    # RNF5 OE
    Candidate_gene_heatmaps(output.dir = DESeq.output.dir, 
                            EnsemblIDs_to_plot = my.candidates, 
                            VST.data = DESeq.RNF5[[1]],
                            my.gene.symbols = my.candidates.symbols, 
                            control_label = 'Empty', 
                            control_n = 4, 
                            control_color = 'black', 
                            treatment_label = 'RNF5_OE', 
                            treatment_n = 4, 
                            treatment_color = 'purple')
    
    
    
    
    
    
    
# Run MDS 

    # IL16
    Run_MDS(VST_expression = DESeq.IL16[[1]], 
            control_color = 'black', 
            control_reps = 4, 
            treatment_label = 'IL16', 
            treatment_color = 'red1', 
            treatment_reps = 4, 
            MDS_dir = MDS.output.dir)
    
    # STARD5
    Run_MDS(VST_expression = DESeq.STARD5[[1]], 
            control_color = 'black', 
            control_reps = 4, 
            treatment_label = 'STARD5', 
            treatment_color = 'pink', 
            treatment_reps = 4, 
            MDS_dir = MDS.output.dir)
    
    # HSD17B12
    Run_MDS(VST_expression = DESeq.HSD17B12[[1]], 
            control_color = 'black', 
            control_reps = 4, 
            treatment_label = 'HSD17B12', 
            treatment_color = '#369193', 
            treatment_reps = 4, 
            MDS_dir = MDS.output.dir)
    
    # RNF5
    Run_MDS(VST_expression = DESeq.RNF5[[1]], 
            control_color = 'black', 
            control_reps = 4, 
            treatment_label = 'RNF5', 
            treatment_color = 'purple', 
            treatment_reps = 4, 
            MDS_dir = MDS.output.dir)



# Run PCA 

    # IL16
    Run_PCA(VST_expression = DESeq.IL16[[1]], 
            control_color = 'black', 
            control_reps = 4, 
            treatment_label = 'IL16', 
            treatment_color = 'red1', 
            treatment_reps = 4, 
            MDS_dir = MDS.output.dir)
    
    # STARD5
    Run_PCA(VST_expression = DESeq.STARD5[[1]], 
            control_color = 'black', 
            control_reps = 4, 
            treatment_label = 'STARD5', 
            treatment_color = 'pink', 
            treatment_reps = 4, 
            MDS_dir = MDS.output.dir)
    
    # HSD17B12
    Run_PCA(VST_expression = DESeq.HSD17B12[[1]], 
            control_color = 'black', 
            control_reps = 4, 
            treatment_label = 'HSD17B12', 
            treatment_color = '#369193', 
            treatment_reps = 4, 
            MDS_dir = MDS.output.dir)
    
    # RNF5
    Run_PCA(VST_expression = DESeq.RNF5[[1]], 
            control_color = 'black', 
            control_reps = 4, 
            treatment_label = 'RNF5', 
            treatment_color = 'purple', 
            treatment_reps = 4, 
            MDS_dir = MDS.output.dir)
    



  
# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/3_DESeq/Session_Info/'
    
    
sink(file = paste(dir.session_info,"Session_Info_GM12878_OE_DESeq.txt", sep =""))
sessionInfo()
sink()      
    



    
# Clean the environment
rm(list=ls())






