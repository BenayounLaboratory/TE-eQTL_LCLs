# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Run_GSEA_Against_Top_eQTLs_Functions.R")

# Load libraries
library(clusterProfiler) # to run GSEA
library(ggplot2) # for split dotplots
library(scales) # for modifying the ggplot colorbar
library(enrichplot) # for gseaplot2








# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/2_GSEA/Results_TElocal_New_eQTLs/'

# Load gene sets 
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Gene_Set_Collections_for_GSEA.R')

# Make gene lists
genelist.EUR.ZSCAN <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/1_DESeq_Genotypes/DESeq_Results/EUR_rs1361387_All_Genes.txt")



# Run GSEA using GO BP Collection
GOBP.ZSCAN <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.EUR.ZSCAN, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'EUR_rs1361387_ZSCAN', output.object = 'filtered')

# Run GSEA using REACTOME Collection
Reactome.ZSCAN <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.EUR.ZSCAN, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'EUR_rs1361387_ZSCAN', output.object = 'filtered')

# Run GSEA using Hallmark Collection
Hallmark.ZSCAN <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.EUR.ZSCAN, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'EUR_rs1361387_ZSCAN', output.object = 'filtered')

# Run GSEA using GTRD TFT Collection
TFT.ZSCAN <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', my.genelist = genelist.EUR.ZSCAN, my.gs.collection = GTRD_Geneset, gs.label = 'GTRD_TFT', condition_label = 'EUR_rs1361387_ZSCAN', output.object = 'filtered')

# Run GSEA using miRDB Collection
miRDB.ZSCAN <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'miRDB/', my.genelist = genelist.EUR.ZSCAN, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = 'EUR_rs1361387_ZSCAN', output.object = 'filtered')



    
# Save results
save(GOBP.ZSCAN, 
     Reactome.ZSCAN, 
     Hallmark.ZSCAN,
     TFT.ZSCAN, 
     miRDB.ZSCAN,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5.R", sep = ''))







# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/2_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_GSEA_TElocal_New_eQTL_Pathways.txt", sep =""))
sessionInfo()
sink()      
    



    
# Clean the environment
rm(list=ls())
    
    
