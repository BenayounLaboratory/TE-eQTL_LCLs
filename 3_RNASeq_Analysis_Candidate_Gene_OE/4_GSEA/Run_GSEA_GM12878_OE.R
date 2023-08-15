# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Run_GSEA_Against_Top_eQTLs_Functions.R") # ***NOTE: The current script relies on functions prepared in the "2_RNASeq_Analysis_Top_eQTLs" analysis folder

# Load libraries
library(clusterProfiler) # to run GSEA
library(ggplot2) # for split dotplots
library(scales) # for modifying the ggplot colorbar
library(enrichplot) # for gseaplot2








# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/4_GSEA/Results_GM12878_Candidate_OE/'

# Load gene sets. ***NOTE: The current script relies on gene sets prepared in the "2_RNASeq_Analysis_Top_eQTLs" analysis folder
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Gene_Set_Collections_for_GSEA.R')

# Make gene lists
genelist.IL16 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/3_DESeq/DESeq_Results_GM12878_OE/All_Genes_IL16.txt")
genelist.STARD5 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/3_DESeq/DESeq_Results_GM12878_OE/All_Genes_STARD5.txt")
genelist.HSD17B12 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/3_DESeq/DESeq_Results_GM12878_OE/All_Genes_HSD17B12.txt")
genelist.RNF5 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/3_DESeq/DESeq_Results_GM12878_OE/All_Genes_RNF5.txt")
  




# Run GSEA using GO BP Collection
GOBP.OE.IL16     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.IL16, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'IL16_OE', output.object = 'filtered')
GOBP.OE.STARD5   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.STARD5, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'STARD5_OE', output.object = 'filtered')
GOBP.OE.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.HSD17B12, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'HSD17B12_OE', output.object = 'filtered')
GOBP.OE.RNF5     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.RNF5, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'RNF5_OE', output.object = 'filtered')

# Run GSEA using REACTOME Collection
Reactome.OE.IL16     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.IL16, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'IL16_OE', output.object = 'filtered')
Reactome.OE.STARD5   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.STARD5, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'STARD5_OE', output.object = 'filtered')
Reactome.OE.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.HSD17B12, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'HSD17B12_OE', output.object = 'filtered')
Reactome.OE.RNF5     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.RNF5, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'RNF5_OE', output.object = 'filtered')

# Run GSEA using Hallmark Collection
Hallmark.OE.IL16     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.IL16, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'IL16_OE', output.object = 'filtered')
Hallmark.OE.STARD5   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.STARD5, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'STARD5_OE', output.object = 'filtered')
Hallmark.OE.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.HSD17B12, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'HSD17B12_OE', output.object = 'filtered')
Hallmark.OE.RNF5     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.RNF5, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'RNF5_OE', output.object = 'filtered')

# Run GSEA using GTRD TFT Collection
TFT.OE.IL16     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', my.genelist = genelist.IL16, my.gs.collection = GTRD_Geneset, gs.label = 'GTRD_TFT', condition_label = 'IL16_OE', output.object = 'filtered')
TFT.OE.STARD5   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', my.genelist = genelist.STARD5, my.gs.collection = GTRD_Geneset, gs.label = 'GTRD_TFT', condition_label = 'STARD5_OE', output.object = 'filtered')
TFT.OE.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', my.genelist = genelist.HSD17B12, my.gs.collection = GTRD_Geneset, gs.label = 'GTRD_TFT', condition_label = 'HSD17B12_OE', output.object = 'filtered')
TFT.OE.RNF5     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', my.genelist = genelist.RNF5, my.gs.collection = GTRD_Geneset, gs.label = 'GTRD_TFT', condition_label = 'RNF5_OE', output.object = 'filtered')

# Run GSEA using miRDB Collection
miRDB.OE.IL16     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'miRDB/', my.genelist = genelist.IL16, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = 'IL16_OE', output.object = 'filtered')
miRDB.OE.STARD5   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'miRDB/', my.genelist = genelist.STARD5, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = 'STARD5_OE', output.object = 'filtered')
miRDB.OE.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'miRDB/', my.genelist = genelist.HSD17B12, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = 'HSD17B12_OE', output.object = 'filtered')
miRDB.OE.RNF5     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'miRDB/', my.genelist = genelist.RNF5, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = 'RNF5_OE', output.object = 'filtered')

# Run GSEA using Repeat Class Collection
Classes.OE.IL16     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.IL16, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'IL16_OE', output.object = 'unfiltered')
Classes.OE.STARD5   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.STARD5, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'STARD5_OE', output.object = 'unfiltered')
Classes.OE.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.HSD17B12, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'HSD17B12_OE', output.object = 'unfiltered')
Classes.OE.RNF5     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.RNF5, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'RNF5_OE', output.object = 'unfiltered')
    
# Run GSEA using Repeat Family Collection
Family.OE.IL16     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.IL16, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'IL16_OE', output.object = 'unfiltered')
Family.OE.STARD5   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.STARD5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'STARD5_OE', output.object = 'unfiltered')
Family.OE.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.HSD17B12, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'HSD17B12_OE', output.object = 'unfiltered')
Family.OE.RNF5     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.RNF5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'RNF5_OE', output.object = 'unfiltered')

    # Generate GSEA Plots for top 3 families

    # IL16
    pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_IL16_OE", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.OE.IL16, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # STARD5
    pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_STARD5_OE", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.OE.STARD5, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # RNF5
    pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_RNF5_OE", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.OE.RNF5, geneSetID = c("Alu subfamilies"), color = c('red2'), base_size = 9, pvalue_table = TRUE)
    dev.off()

    
    
    
# Run GSEA using L1 Subfamily Collection
L1.OE.IL16     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_L1s_Only/', my.genelist = genelist.IL16, my.gs.collection = L1_by_age.gs, gs.label = 'L1_Subfamilies', condition_label = 'IL16_OE', output.object = 'unfiltered')
L1.OE.STARD5   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_L1s_Only/', my.genelist = genelist.STARD5, my.gs.collection = L1_by_age.gs, gs.label = 'L1_Subfamilies', condition_label = 'STARD5_OE', output.object = 'unfiltered')
L1.OE.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_L1s_Only/', my.genelist = genelist.HSD17B12, my.gs.collection = L1_by_age.gs, gs.label = 'L1_Subfamilies', condition_label = 'HSD17B12_OE', output.object = 'unfiltered')
L1.OE.RNF5     <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_L1s_Only/', my.genelist = genelist.RNF5, my.gs.collection = L1_by_age.gs, gs.label = 'L1_Subfamilies', condition_label = 'RNF5_OE', output.object = 'unfiltered')

    # Generate GSEA Plots

    # IL16
    pdf(paste(my.output.dir, "Repeat_L1s_Only/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_L1_Subfamilies_IL16_OE", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(L1.OE.IL16, geneSetID = c("L1M subfamilies", "L1P subfamilies", "L1PA subfamilies"), color = c('coral', 'firebrick2', 'firebrick4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # STARD5
    pdf(paste(my.output.dir, "Repeat_L1s_Only/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_L1_Subfamilies_STARD5_OE", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(L1.OE.STARD5, geneSetID = c("L1M subfamilies", "L1P subfamilies", "L1PA subfamilies"), color = c('coral', 'firebrick2', 'firebrick4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # RNF5
    pdf(paste(my.output.dir, "Repeat_L1s_Only/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_L1_Subfamilies_RNF5_OE", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(L1.OE.RNF5, geneSetID = c("L1M subfamilies", "L1P subfamilies", "L1PA subfamilies"), color = c('coral', 'firebrick2', 'firebrick4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    
    
    



    
# Save results
save(GOBP.OE.IL16, GOBP.OE.STARD5, GOBP.OE.HSD17B12, GOBP.OE.RNF5,
     Reactome.OE.IL16, Reactome.OE.STARD5, Reactome.OE.HSD17B12, Reactome.OE.RNF5,
     Hallmark.OE.IL16, Hallmark.OE.STARD5, Hallmark.OE.HSD17B12, Hallmark.OE.RNF5,
     TFT.OE.IL16, TFT.OE.STARD5, TFT.OE.HSD17B12, TFT.OE.RNF5,
     miRDB.OE.IL16, miRDB.OE.STARD5, miRDB.OE.HSD17B12, miRDB.OE.RNF5,
     Classes.OE.IL16, Classes.OE.STARD5, Classes.OE.HSD17B12, Classes.OE.RNF5,
     Family.OE.IL16, Family.OE.STARD5, Family.OE.HSD17B12, Family.OE.RNF5,
     L1.OE.IL16, L1.OE.STARD5, L1.OE.HSD17B12, L1.OE.RNF5,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5.R", sep = ''))







# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/4_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_GM12878_OE_GSEA.txt", sep =""))
sessionInfo()
sink()      
    



    
# Clean the environment
rm(list=ls())
    
    


