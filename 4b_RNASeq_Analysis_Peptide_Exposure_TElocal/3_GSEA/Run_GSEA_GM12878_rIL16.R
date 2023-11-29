# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Run_GSEA_Against_Top_eQTLs_Functions.R")

# Load libraries
library(clusterProfiler) # to run GSEA
library(ggplot2) # for split dotplots
library(scales) # for modifying the ggplot colorbar
library(enrichplot) # for gseaplot2






# GENERAL PARAMETERS

# Load gene sets 
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/2_GSEA/Repeat_Subset_Gene_Set_Collections_for_GSEA.R')

# Make gene lists
genelist.rIL16.T24 <-   make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/2_DESeq/DESeq_Results_GM12878_rIL16/rIL16_T24_vs_T0_All_Genes.txt")
genelist.rIL16.T48 <-   make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/2_DESeq/DESeq_Results_GM12878_rIL16/rIL16_T48_vs_T0_All_Genes.txt")





# RUN GSEA USING INTRONIC TE GENE SETS ------------------------------------------------------------------------------------------------------------

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/3_GSEA/Results_GM12878_rIL16/Repeat_Families_intronic/'

# Define the gene set of interest
Repeat_family.gs <- intronic_Repeat_family.gs

# Run GSEA using Repeat Family Collection
intronic.T24   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.rIL16.T24, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'rIL16_T24', output.object = 'unfiltered')
intronic.T48   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.rIL16.T48, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'rIL16_T48', output.object = 'unfiltered')

    # # Generate GSEA Plots for top 3 families
    # 
    # # T24
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_rIL16_T24", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.T24, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # T48
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_rIL16_T48", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.T48, geneSetID = c("L1 subfamilies"), color = c( 'red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    
    
    
# Save results
save(intronic.T24, intronic.T48,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5_Intronic_Repeat_Families.R", sep = ''))   



# RUN GSEA USING INTERGENIC DISTAL TE GENE SETS ------------------------------------------------------------------------------------------------------------

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/3_GSEA/Results_GM12878_rIL16/Repeat_Families_intergenic_distal/'

# Define the gene set of interest
Repeat_family.gs <- intergenic_distal_Repeat_family.gs

# Run GSEA using Repeat Family Collection
distal.T24   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.rIL16.T24, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'rIL16_T24', output.object = 'unfiltered')
distal.T48   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.rIL16.T48, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'rIL16_T48', output.object = 'unfiltered')

    # # Generate GSEA Plots for top 3 families
    # 
    # # T24
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_rIL16_T24", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.T24, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # T48
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_rIL16_T48", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.T48, geneSetID = c("L1 subfamilies"), color = c( 'red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    
    
    
# Save results
save(distal.T24, distal.T48,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5_Intergenic_Distal_Repeat_Families.R", sep = ''))   




# RUN GSEA USING INTERGENIC NEARBY TE GENE SETS ------------------------------------------------------------------------------------------------------------

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/3_GSEA/Results_GM12878_rIL16/Repeat_Families_intergenic_nearby/'

# Define the gene set of interest
Repeat_family.gs <- intergenic_nearby_Repeat_family.gs

# Run GSEA using Repeat Family Collection
nearby.T24   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.rIL16.T24, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'rIL16_T24', output.object = 'unfiltered')
nearby.T48   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.rIL16.T48, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'rIL16_T48', output.object = 'unfiltered')

    # # Generate GSEA Plots for top 3 families
    # 
    # # T24
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_rIL16_T24", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.T24, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # T48
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_rIL16_T48", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.T48, geneSetID = c("L1 subfamilies"), color = c( 'red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    
    
    
# Save results
save(nearby.T24, nearby.T48,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5_Intergenic_Nearby_Repeat_Families.R", sep = ''))   





# RUN GSEA USING EXONIC TE GENE SETS ------------------------------------------------------------------------------------------------------------

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/3_GSEA/Results_GM12878_rIL16/Repeat_Families_exonic/'

# Define the gene set of interest
Repeat_family.gs <- exonic_Repeat_family.gs

# Run GSEA using Repeat Family Collection
exonic.T24   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.rIL16.T24, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'rIL16_T24', output.object = 'unfiltered')
exonic.T48   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.rIL16.T48, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'rIL16_T48', output.object = 'unfiltered')

    # # Generate GSEA Plots for top 3 families
    # 
    # # T24
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_rIL16_T24", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.T24, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # T48
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_rIL16_T48", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.T48, geneSetID = c("L1 subfamilies"), color = c( 'red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    
    
    
# Save results
save(exonic.T24, exonic.T48,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5_Exonic_Repeat_Families.R", sep = ''))   





# SESSION INFO ------------------------------------------------------------------------------------------------------------
 
    
# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/3_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_GSEA_GM12878_rIL16.txt", sep =""))
sessionInfo()
sink()      
    



    
# Clean the environment
rm(list=ls())
    
    


