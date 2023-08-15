# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Run_GSEA_Against_Top_eQTLs_Functions.R")

# Load libraries
library(clusterProfiler) # to run GSEA
library(ggplot2) # for split dotplots
library(scales) # for modifying the ggplot colorbar
library(enrichplot) # for gseaplot2





# RUN GSEA ----------------------------------------------------------------------



# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/4_GSEA/Results_GM12878_rIL16/'

# Make gene lists
genelist.rIL16.T24 <-   make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/3_DESeq/DESeq_Results_GM12878_rIL16/rIL16_T24_vs_T0_All_Genes.txt")
genelist.rIL16.T48 <-   make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/3_DESeq/DESeq_Results_GM12878_rIL16/rIL16_T48_vs_T0_All_Genes.txt")
genelist.rIL16.delta <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/3_DESeq/DESeq_Results_GM12878_rIL16/rIL16_T48_vs_T24_All_Genes.txt")

# Load gene sets ***NOTE: The current script relies on gene sets prepared in the "2_RNASeq_Analysis_Top_eQTLs" analysis folder
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Gene_Set_Collections_for_GSEA.R')



# Run GSEA using GO BP Collection
GOBP.T24   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.rIL16.T24, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'rIL16_T24', output.object = 'filtered')
GOBP.T48   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.rIL16.T48, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'rIL16_T48', output.object = 'filtered')
GOBP.delta <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.rIL16.delta, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'rIL16_delta', output.object = 'filtered')

# Run GSEA using REACTOME Collection
Reactome.T24   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.rIL16.T24, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'rIL16_T24', output.object = 'filtered')
Reactome.T48   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.rIL16.T48, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'rIL16_T48', output.object = 'filtered')
Reactome.delta <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.rIL16.delta, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'rIL16_delta', output.object = 'filtered')

# Run GSEA using Hallmark Collection
Hallmark.T24   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.rIL16.T24, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'rIL16_T24', output.object = 'filtered')
Hallmark.T48   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.rIL16.T48, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'rIL16_T48', output.object = 'filtered')
Hallmark.delta <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.rIL16.delta, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'rIL16_delta', output.object = 'filtered')

# Run GSEA using GTRD TFT Collection
TFT.T24   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', my.genelist = genelist.rIL16.T24, my.gs.collection = GTRD_Geneset, gs.label = 'GTRD_TFT', condition_label = 'rIL16_T24', output.object = 'filtered')
TFT.T48   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', my.genelist = genelist.rIL16.T48, my.gs.collection = GTRD_Geneset, gs.label = 'GTRD_TFT', condition_label = 'rIL16_T48', output.object = 'filtered')
TFT.delta <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', my.genelist = genelist.rIL16.delta, my.gs.collection = GTRD_Geneset, gs.label = 'GTRD_TFT', condition_label = 'rIL16_delta', output.object = 'filtered')

# Run GSEA using miRDB Collection
miRNA.T24   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'miRDB/', my.genelist = genelist.rIL16.T24, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = 'rIL16_T24', output.object = 'filtered')
miRNA.T48   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'miRDB/', my.genelist = genelist.rIL16.T48, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = 'rIL16_T48', output.object = 'filtered')
miRNA.delta <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'miRDB/', my.genelist = genelist.rIL16.delta, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = 'rIL16_delta', output.object = 'filtered')

    # Generate GSEA plots for known miRNA regulators (let7, mir128 )

    # T24
    pdf(paste(my.output.dir, "miRDB/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_let7_rIL16_T24", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(miRNA.T24, geneSetID = c("LET-7F-2-3P"), title = "let-7 Targets GSEA", color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # delta
    pdf(paste(my.output.dir, "miRDB/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_let7_rIL16_delta", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(miRNA.delta, geneSetID = c("LET-7F-2-3P"), title = "let-7 Targets GSEA", color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()

# Run GSEA using Repeat Class Collection
Classes.T24   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.rIL16.T24, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'rIL16_T24', output.object = 'unfiltered')
Classes.T48   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.rIL16.T48, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'rIL16_T48', output.object = 'unfiltered')
Classes.delta <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.rIL16.delta, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'rIL16_delta', output.object = 'unfiltered')
    
# Run GSEA using Repeat Family Collection
Family.T24   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.rIL16.T24, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'rIL16_T24', output.object = 'unfiltered')
Family.T48   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.rIL16.T48, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'rIL16_T48', output.object = 'unfiltered')
Family.delta <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.rIL16.delta, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'rIL16_delta', output.object = 'unfiltered')

    # Generate GSEA Plots for top 3 families

    # T24
    pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_rIL16_T24", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.T24, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # T48
    pdf(paste(my.output.dir, "Repeat_Families/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_Repeat_Families_rIL16_T48", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.T48, geneSetID = c("L1 subfamilies"), color = c( 'red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    
# Run GSEA using L1 Age-segregated Collection
L1.T24   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_L1s_Only/', my.genelist = genelist.rIL16.T24, my.gs.collection = L1_by_age.gs, gs.label = 'L1_Subfamilies', condition_label = 'rIL16_T24', output.object = 'unfiltered')
L1.T48   <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_L1s_Only/', my.genelist = genelist.rIL16.T48, my.gs.collection = L1_by_age.gs, gs.label = 'L1_Subfamilies', condition_label = 'rIL16_T48', output.object = 'unfiltered')
L1.delta <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_L1s_Only/', my.genelist = genelist.rIL16.delta, my.gs.collection = L1_by_age.gs, gs.label = 'L1_Subfamilies', condition_label = 'rIL16_delta', output.object = 'unfiltered')


    # Generate GSEA Plots 

    # T24
    pdf(paste(my.output.dir, "Repeat_L1s_Only/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_L1_Subfamilies_rIL16_T24", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(L1.T24, geneSetID = c("L1M subfamilies", "L1P subfamilies", "L1PA subfamilies"), color = c('coral', 'firebrick2', 'firebrick4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # T48
    pdf(paste(my.output.dir, "Repeat_L1s_Only/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_L1_Subfamilies_rIL16_T48", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(L1.T48, geneSetID = c("L1M subfamilies", "L1P subfamilies", "L1PA subfamilies"), color = c('coral', 'firebrick2', 'firebrick4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # delta
    pdf(paste(my.output.dir, "Repeat_L1s_Only/", "ClusterProfiler_GSEA_EnrichmentPlot_FDR5_L1_Subfamilies_rIL16_delta", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(L1.delta, geneSetID = c("L1M subfamilies", "L1P subfamilies", "L1PA subfamilies"), color = c('coral', 'firebrick2', 'firebrick4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    
    
    



    
# Save results
save(GOBP.T24, GOBP.T48, GOBP.delta,
     Reactome.T24, Reactome.T48, Reactome.delta,
     Hallmark.T24, Hallmark.T48, Hallmark.delta,
     TFT.T24, TFT.T48, TFT.delta,
     miRNA.T24, miRNA.T48, miRNA.delta,
     Classes.T24, Classes.T48, Classes.delta,
     Family.T24, Family.T48, Family.delta,
     L1.T24, L1.T48, L1.delta,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5.R", sep = ''))    
    



    
    
# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/4_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_GSEA_GM12878_rIL16.txt", sep =""))
sessionInfo()
sink()      
    



    
# Clean the environment
rm(list=ls())
    
    


