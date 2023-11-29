# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Run_GSEA_Against_Top_eQTLs_Functions.R") # ***NOTE: The current script relies on functions prepared in the "2_RNASeq_Analysis_Top_eQTLs" analysis folder

# Load libraries
library(clusterProfiler) # to run GSEA
library(ggplot2) # for split dotplots
library(scales) # for modifying the ggplot colorbar
library(enrichplot) # for gseaplot2






# GENERAL PARAMETERS

# Load gene sets 
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/2_GSEA/Repeat_Subset_Gene_Set_Collections_for_GSEA.R')

# Make gene lists
genelist.IL16 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3b_RNASeq_Analysis_Candidate_Gene_OE_TElocal/2_DESeq/DESeq_Results_GM12878_OE/All_Genes_IL16.txt")
genelist.STARD5 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3b_RNASeq_Analysis_Candidate_Gene_OE_TElocal/2_DESeq/DESeq_Results_GM12878_OE/All_Genes_STARD5.txt")
genelist.HSD17B12 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3b_RNASeq_Analysis_Candidate_Gene_OE_TElocal/2_DESeq/DESeq_Results_GM12878_OE/All_Genes_HSD17B12.txt")
genelist.RNF5 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3b_RNASeq_Analysis_Candidate_Gene_OE_TElocal/2_DESeq/DESeq_Results_GM12878_OE/All_Genes_RNF5.txt")
  


# RUN GSEA USING INTRONIC TE GENE SETS ------------------------------------------------------------------------------------------------------------

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3b_RNASeq_Analysis_Candidate_Gene_OE_TElocal/3_GSEA/Results_GM12878_Candidate_OE/Repeat_Families_intronic/'

# Define the gene set of interest
Repeat_family.gs <- intronic_Repeat_family.gs

# Run GSEA using Repeat Family Collection
intronic.OE.IL16     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.IL16, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'IL16_OE', output.object = 'unfiltered')
intronic.OE.STARD5   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.STARD5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'STARD5_OE', output.object = 'unfiltered')
intronic.OE.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.HSD17B12, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'HSD17B12_OE', output.object = 'unfiltered')
intronic.OE.RNF5     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.RNF5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'RNF5_OE', output.object = 'unfiltered')

    # # Generate GSEA Plots for top 3 families
    # 
    # # IL16
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_IL16_OE", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.OE.IL16, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # STARD5
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_STARD5_OE", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.OE.STARD5, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # RNF5
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_RNF5_OE", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.OE.RNF5, geneSetID = c("Alu subfamilies"), color = c('red2'), base_size = 9, pvalue_table = TRUE)
    # dev.off()

    
# Save results
save(intronic.OE.IL16, intronic.OE.STARD5, intronic.OE.HSD17B12, intronic.OE.RNF5,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5_Intronic_Repeat_Families.R", sep = ''))    
    

# RUN GSEA USING INTERGENIC DISTAL TE GENE SETS ------------------------------------------------------------------------------------------------------------

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3b_RNASeq_Analysis_Candidate_Gene_OE_TElocal/3_GSEA/Results_GM12878_Candidate_OE/Repeat_Families_intergenic_distal/'

# Define the gene set of interest
Repeat_family.gs <- intergenic_distal_Repeat_family.gs

# Run GSEA using Repeat Family Collection
distal.OE.IL16     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.IL16, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'IL16_OE', output.object = 'unfiltered')
distal.OE.STARD5   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.STARD5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'STARD5_OE', output.object = 'unfiltered')
distal.OE.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.HSD17B12, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'HSD17B12_OE', output.object = 'unfiltered')
distal.OE.RNF5     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.RNF5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'RNF5_OE', output.object = 'unfiltered')

    # # Generate GSEA Plots for top 3 families
    # 
    # # IL16
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_IL16_OE", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.OE.IL16, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # STARD5
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_STARD5_OE", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.OE.STARD5, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # RNF5
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_RNF5_OE", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.OE.RNF5, geneSetID = c("Alu subfamilies"), color = c('red2'), base_size = 9, pvalue_table = TRUE)
    # dev.off()

    
# Save results
save(distal.OE.IL16, distal.OE.STARD5, distal.OE.HSD17B12, distal.OE.RNF5,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5_Intergenic_Distal_Repeat_Families.R", sep = ''))    
    


# RUN GSEA USING INTERGENIC NEARBY TE GENE SETS ------------------------------------------------------------------------------------------------------------

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3b_RNASeq_Analysis_Candidate_Gene_OE_TElocal/3_GSEA/Results_GM12878_Candidate_OE/Repeat_Families_intergenic_nearby/'

# Define the gene set of interest
Repeat_family.gs <- intergenic_nearby_Repeat_family.gs

# Run GSEA using Repeat Family Collection
nearby.OE.IL16     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.IL16, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'IL16_OE', output.object = 'unfiltered')
nearby.OE.STARD5   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.STARD5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'STARD5_OE', output.object = 'unfiltered')
nearby.OE.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.HSD17B12, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'HSD17B12_OE', output.object = 'unfiltered')
nearby.OE.RNF5     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.RNF5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'RNF5_OE', output.object = 'unfiltered')

    # # Generate GSEA Plots for top 3 families
    # 
    # # IL16
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_IL16_OE", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.OE.IL16, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # STARD5
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_STARD5_OE", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.OE.STARD5, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # RNF5
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_RNF5_OE", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.OE.RNF5, geneSetID = c("Alu subfamilies"), color = c('red2'), base_size = 9, pvalue_table = TRUE)
    # dev.off()

    
# Save results
save(nearby.OE.IL16, nearby.OE.STARD5, nearby.OE.HSD17B12, nearby.OE.RNF5,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5_Intergenic_Nearby_Repeat_Families.R", sep = ''))    
    



# RUN GSEA USING EXONIC TE GENE SETS ------------------------------------------------------------------------------------------------------------

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3b_RNASeq_Analysis_Candidate_Gene_OE_TElocal/3_GSEA/Results_GM12878_Candidate_OE/Repeat_Families_exonic/'

# Define the gene set of interest
Repeat_family.gs <- exonic_Repeat_family.gs

# Run GSEA using Repeat Family Collection
exonic.OE.IL16     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.IL16, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'IL16_OE', output.object = 'unfiltered')
exonic.OE.STARD5   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.STARD5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'STARD5_OE', output.object = 'unfiltered')
exonic.OE.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.HSD17B12, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'HSD17B12_OE', output.object = 'unfiltered')
exonic.OE.RNF5     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.RNF5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'RNF5_OE', output.object = 'unfiltered')

    # # Generate GSEA Plots for top 3 families
    # 
    # # IL16
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_IL16_OE", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.OE.IL16, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # STARD5
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_STARD5_OE", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.OE.STARD5, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # RNF5
    # pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_RNF5_OE", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.OE.RNF5, geneSetID = c("Alu subfamilies"), color = c('red2'), base_size = 9, pvalue_table = TRUE)
    # dev.off()

    
# Save results
save(exonic.OE.IL16, exonic.OE.STARD5, exonic.OE.HSD17B12, exonic.OE.RNF5,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5_Exonic_Repeat_Families.R", sep = ''))    
    




# SESSION INFO ------------------------------------------------------------------------------------------------------------


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3b_RNASeq_Analysis_Candidate_Gene_OE_TElocal/3_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_GM12878_OE_GSEA.txt", sep =""))
sessionInfo()
sink()      
    



    
# Clean the environment
rm(list=ls())
    
    


