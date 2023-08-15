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
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Results_Against_Top_eQTLs/'

# Load gene sets 
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Gene_Set_Collections_for_GSEA.R')

# Make gene lists
genelist.EUR.IL16 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/1_DESeq_Genotypes/DESeq_Results/EUR_rs11635336_All_Genes.txt")
genelist.EUR.HLA <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/1_DESeq_Genotypes/DESeq_Results/EUR_rs9271894_All_Genes.txt")
genelist.EUR.HSD17B12 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/1_DESeq_Genotypes/DESeq_Results/EUR_rs1061810_All_Genes.txt")
genelist.EUR.SNV4 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/1_DESeq_Genotypes/DESeq_Results/EUR_rs112581165_All_Genes.txt")
genelist.EUR.SNV5 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/1_DESeq_Genotypes/DESeq_Results/EUR_rs72691418_All_Genes.txt")
genelist.EUR.RNF5 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/1_DESeq_Genotypes/DESeq_Results/EUR_rs9270493_All_Genes.txt")


genelist.YRI.HSD17B12 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/1_DESeq_Genotypes/DESeq_Results/YRI_rs2176598_All_Genes.txt")
genelist.YRI.HLA <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/1_DESeq_Genotypes/DESeq_Results/YRI_rs9271379_All_Genes.txt")
  






# Run GSEA using GO BP Collection
GOBP.IL16 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.EUR.IL16, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'EUR_rs11635336_IL16', output.object = 'filtered')
GOBP.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.EUR.HLA, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'EUR_rs9271894_HLA', output.object = 'filtered')
GOBP.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.EUR.HSD17B12, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'EUR_rs1061810_HSD17B12', output.object = 'filtered')
GOBP.SNV4 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.EUR.SNV4, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'EUR_rs112581165_SNV4', output.object = 'filtered')
GOBP.SNV5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.EUR.SNV5, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'EUR_rs72691418_SNV5', output.object = 'filtered')
GOBP.RNF5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.EUR.RNF5, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'EUR_rs9270493_RNF5', output.object = 'filtered')

GOBP.YRI.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.YRI.HSD17B12, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'YRI_rs2176598_HSD17B12', output.object = 'filtered')
GOBP.YRI.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'GO/', my.genelist = genelist.YRI.HLA, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = 'YRI_rs9271379_HLA', output.object = 'filtered')

# Run GSEA using REACTOME Collection
Reactome.IL16 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.EUR.IL16, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'EUR_rs11635336_IL16', output.object = 'filtered')
Reactome.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.EUR.HLA, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'EUR_rs9271894_HLA', output.object = 'filtered')
Reactome.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.EUR.HSD17B12, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'EUR_rs1061810_HSD17B12', output.object = 'filtered')
Reactome.SNV4 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.EUR.SNV4, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'EUR_rs112581165_SNV4', output.object = 'filtered')
Reactome.SNV5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.EUR.SNV5, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'EUR_rs72691418_SNV5', output.object = 'filtered')
Reactome.RNF5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.EUR.RNF5, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'EUR_rs9270493_RNF5', output.object = 'filtered')

Reactome.YRI.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.YRI.HSD17B12, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'YRI_rs2176598_HSD17B12', output.object = 'filtered')
Reactome.YRI.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Reactome/', my.genelist = genelist.YRI.HLA, my.gs.collection = Reactome_Geneset, gs.label = 'REACTOME', condition_label = 'YRI_rs9271379_HLA', output.object = 'filtered')

# Run GSEA using Hallmark Collection
Hallmark.IL16 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.EUR.IL16, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'EUR_rs11635336_IL16', output.object = 'filtered')
Hallmark.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.EUR.HLA, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'EUR_rs9271894_HLA', output.object = 'filtered')
Hallmark.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.EUR.HSD17B12, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'EUR_rs1061810_HSD17B12', output.object = 'filtered')
Hallmark.SNV4 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.EUR.SNV4, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'EUR_rs112581165_SNV4', output.object = 'filtered')
Hallmark.SNV5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.EUR.SNV5, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'EUR_rs72691418_SNV5', output.object = 'filtered')
Hallmark.RNF5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.EUR.RNF5, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'EUR_rs9270493_RNF5', output.object = 'filtered')

Hallmark.YRI.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.YRI.HSD17B12, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'YRI_rs2176598_HSD17B12', output.object = 'filtered')
Hallmark.YRI.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', my.genelist = genelist.YRI.HLA, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = 'YRI_rs9271379_HLA', output.object = 'filtered')

# Run GSEA using GTRD TFT Collection
TFT.IL16 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', my.genelist = genelist.EUR.IL16, my.gs.collection = GTRD_Geneset, gs.label = 'GTRD_TFT', condition_label = 'EUR_rs11635336_IL16', output.object = 'filtered')
TFT.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', my.genelist = genelist.EUR.HLA, my.gs.collection = GTRD_Geneset, gs.label = 'GTRD_TFT', condition_label = 'EUR_rs9271894_HLA', output.object = 'filtered')
TFT.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', my.genelist = genelist.EUR.HSD17B12, my.gs.collection = GTRD_Geneset, gs.label = 'GTRD_TFT', condition_label = 'EUR_rs1061810_HSD17B12', output.object = 'filtered')
TFT.SNV4 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', my.genelist = genelist.EUR.SNV4, my.gs.collection = GTRD_Geneset, gs.label = 'GTRD_TFT', condition_label = 'EUR_rs112581165_SNV4', output.object = 'filtered')
TFT.SNV5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', my.genelist = genelist.EUR.SNV5, my.gs.collection = GTRD_Geneset, gs.label = 'GTRD_TFT', condition_label = 'EUR_rs72691418_SNV5', output.object = 'filtered')
TFT.RNF5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', my.genelist = genelist.EUR.RNF5, my.gs.collection = GTRD_Geneset, gs.label = 'GTRD_TFT', condition_label = 'EUR_rs9270493_RNF5', output.object = 'filtered')

TFT.YRI.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', my.genelist = genelist.YRI.HSD17B12, my.gs.collection = GTRD_Geneset, gs.label = 'GTRD_TFT', condition_label = 'YRI_rs2176598_HSD17B12', output.object = 'filtered')
TFT.YRI.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', my.genelist = genelist.YRI.HLA, my.gs.collection = GTRD_Geneset, gs.label = 'GTRD_TFT', condition_label = 'YRI_rs9271379_HLA', output.object = 'filtered')

# Run GSEA using miRDB Collection
miRDB.IL16 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'miRDB/', my.genelist = genelist.EUR.IL16, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = 'EUR_rs11635336_IL16', output.object = 'filtered')
miRDB.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'miRDB/', my.genelist = genelist.EUR.HLA, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = 'EUR_rs9271894_HLA', output.object = 'filtered')
miRDB.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'miRDB/', my.genelist = genelist.EUR.HSD17B12, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = 'EUR_rs1061810_HSD17B12', output.object = 'filtered')
miRDB.SNV4 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'miRDB/', my.genelist = genelist.EUR.SNV4, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = 'EUR_rs112581165_SNV4', output.object = 'filtered')
miRDB.SNV5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'miRDB/', my.genelist = genelist.EUR.SNV5, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = 'EUR_rs72691418_SNV5', output.object = 'filtered')
miRDB.RNF5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'miRDB/', my.genelist = genelist.EUR.RNF5, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = 'EUR_rs9270493_RNF5', output.object = 'filtered')

miRDB.YRI.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'miRDB/', my.genelist = genelist.YRI.HSD17B12, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = 'YRI_rs2176598_HSD17B12', output.object = 'filtered')
miRDB.YRI.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'miRDB/', my.genelist = genelist.YRI.HLA, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = 'YRI_rs9271379_HLA', output.object = 'filtered')

# Run GSEA using Repeat Class Collection
Class.EUR.IL16 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.EUR.IL16, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'EUR_rs11635336_IL16', output.object = 'unfiltered')
Class.EUR.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.EUR.HLA, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'EUR_rs9271894_HLA', output.object = 'unfiltered')
Class.EUR.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.EUR.HSD17B12, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'EUR_rs1061810_HSD17B12', output.object = 'unfiltered')
Class.EUR.SNV4 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.EUR.SNV4, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'EUR_rs112581165_SNV4', output.object = 'unfiltered')
Class.EUR.SNV5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.EUR.SNV5, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'EUR_rs72691418_SNV5', output.object = 'unfiltered')
Class.EUR.RNF5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.EUR.RNF5, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'EUR_rs9270493_RNF5', output.object = 'unfiltered')

Class.YRI.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.YRI.HSD17B12, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'YRI_rs2176598_HSD17B12', output.object = 'unfiltered')
Class.YRI.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Classes/', my.genelist = genelist.YRI.HLA, my.gs.collection = Repeat_class.gs, gs.label = 'Repeat_Classes', condition_label = 'YRI_rs9271379_HLA', output.object = 'unfiltered') 
    
# Run GSEA using Repeat Family Collection
Family.EUR.IL16 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.EUR.IL16, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs11635336_IL16', output.object = 'unfiltered')
Family.EUR.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.EUR.HLA, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs9271894_HLA', output.object = 'unfiltered')
Family.EUR.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.EUR.HSD17B12, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs1061810_HSD17B12', output.object = 'unfiltered')
Family.EUR.SNV4 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.EUR.SNV4, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs112581165_SNV4', output.object = 'unfiltered')
Family.EUR.SNV5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.EUR.SNV5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs72691418_SNV5', output.object = 'unfiltered')
Family.EUR.RNF5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.EUR.RNF5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs9270493_RNF5', output.object = 'unfiltered')

Family.YRI.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.YRI.HSD17B12, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'YRI_rs2176598_HSD17B12', output.object = 'unfiltered')
Family.YRI.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_Families/', my.genelist = genelist.YRI.HLA, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'YRI_rs9271379_HLA', output.object = 'unfiltered')

    # Generate GSEA Plots for L1 family

    # EUR IL16 SNP
    pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs11635336_IL16", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.EUR.IL16, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # EUR HLA SNP
    pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs9271894_HLA", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.EUR.HLA, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # EUR HSD17B12 SNP
    pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs1061810_HSD17B12", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.EUR.HSD17B12, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # EUR SNV4
    pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs112581165_SNV4", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.EUR.SNV4, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # EUR SNV5
    pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs72691418_SNV5", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.EUR.SNV5, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # EUR RNF5 SNP
    pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs9270493_RNF5", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.EUR.RNF5, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # YRI HSD17B12 SNP
    pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_YRI_rs2176598_HSD17B12", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.YRI.HSD17B12, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()

    # YRI HLA SNP
    pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_YRI_rs9271379_HLA", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(Family.YRI.HLA, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()

    
    
    
# Run GSEA using L1 Subfamily Collection
L1.EUR.IL16 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_L1s_Only/', my.genelist = genelist.EUR.IL16, my.gs.collection = L1_by_age.gs, gs.label = 'L1_Subfamilies', condition_label = 'EUR_rs11635336_IL16', output.object = 'unfiltered')
L1.EUR.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_L1s_Only/', my.genelist = genelist.EUR.HLA, my.gs.collection = L1_by_age.gs, gs.label = 'L1_Subfamilies', condition_label = 'EUR_rs9271894_HLA', output.object = 'unfiltered')
L1.EUR.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_L1s_Only/', my.genelist = genelist.EUR.HSD17B12, my.gs.collection = L1_by_age.gs, gs.label = 'L1_Subfamilies', condition_label = 'EUR_rs1061810_HSD17B12', output.object = 'unfiltered')
L1.EUR.SNV4 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_L1s_Only/', my.genelist = genelist.EUR.SNV4, my.gs.collection = L1_by_age.gs, gs.label = 'L1_Subfamilies', condition_label = 'EUR_rs112581165_SNV4', output.object = 'unfiltered')
L1.EUR.SNV5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_L1s_Only/', my.genelist = genelist.EUR.SNV5, my.gs.collection = L1_by_age.gs, gs.label = 'L1_Subfamilies', condition_label = 'EUR_rs72691418_SNV5', output.object = 'unfiltered')
L1.EUR.RNF5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_L1s_Only/', my.genelist = genelist.EUR.RNF5, my.gs.collection = L1_by_age.gs, gs.label = 'L1_Subfamilies', condition_label = 'EUR_rs9270493_RNF5', output.object = 'unfiltered')

L1.YRI.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_L1s_Only/', my.genelist = genelist.YRI.HSD17B12, my.gs.collection = L1_by_age.gs, gs.label = 'L1_Subfamilies', condition_label = 'YRI_rs2176598_HSD17B12', output.object = 'unfiltered')
L1.YRI.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = 'Repeat_L1s_Only/', my.genelist = genelist.YRI.HLA, my.gs.collection = L1_by_age.gs, gs.label = 'L1_Subfamilies', condition_label = 'YRI_rs9271379_HLA', output.object = 'unfiltered')

    # Generate GSEA Plots for age-segragated L1 Subfamilies

    # EUR IL16 SNP
    pdf(paste(my.output.dir, "Repeat_L1s_Only/", "GSEA_EnrichmentPlot_FDR5_L1_Subamilies_EUR_rs11635336_IL16", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(L1.EUR.IL16, geneSetID = c("L1M subfamilies", "L1P subfamilies", "L1PA subfamilies"), color = c('coral', 'red2', 'red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # EUR HLA SNP
    pdf(paste(my.output.dir, "Repeat_L1s_Only/", "GSEA_EnrichmentPlot_FDR5_L1_Subamilies_EUR_rs9271894_HLA", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(L1.EUR.HLA, geneSetID = c("L1M subfamilies", "L1P subfamilies", "L1PA subfamilies"), color = c('coral', 'red2', 'red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # EUR HSD17B12 SNP
    pdf(paste(my.output.dir, "Repeat_L1s_Only/", "GSEA_EnrichmentPlot_FDR5_L1_Subamilies_EUR_rs1061810_HSD17B12", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(L1.EUR.HSD17B12, geneSetID = c("L1M subfamilies", "L1P subfamilies", "L1PA subfamilies"), color = c('coral', 'red2', 'red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # EUR SNV4
    pdf(paste(my.output.dir, "Repeat_L1s_Only/", "GSEA_EnrichmentPlot_FDR5_L1_Subamilies_EUR_rs112581165_SNV4", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(L1.EUR.SNV4, geneSetID = c("L1M subfamilies", "L1P subfamilies", "L1PA subfamilies"), color = c('coral', 'red2', 'red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # EUR SNV5
    pdf(paste(my.output.dir, "Repeat_L1s_Only/", "GSEA_EnrichmentPlot_FDR5_L1_Subamilies_EUR_rs72691418_SNV5", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(L1.EUR.SNV5, geneSetID = c("L1M subfamilies", "L1P subfamilies", "L1PA subfamilies"), color = c('coral', 'red2', 'red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # EUR RNF5 SNP
    pdf(paste(my.output.dir, "Repeat_L1s_Only/", "GSEA_EnrichmentPlot_FDR5_L1_Subamilies_EUR_rs9270493_RNF5", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(L1.EUR.RNF5, geneSetID = c("L1M subfamilies", "L1P subfamilies", "L1PA subfamilies"), color = c('coral', 'red2', 'red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    # YRI HSD17B12 SNP
    pdf(paste(my.output.dir, "Repeat_L1s_Only/", "GSEA_EnrichmentPlot_FDR5_L1_Subamilies_YRI_rs2176598_HSD17B12", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(L1.YRI.HSD17B12, geneSetID = c("L1M subfamilies", "L1P subfamilies", "L1PA subfamilies"), color = c('coral', 'red2', 'red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()

    # YRI HLA SNP
    pdf(paste(my.output.dir, "Repeat_L1s_Only/", "GSEA_EnrichmentPlot_L1_Subamilies_YRI_rs9271379_HLA", ".pdf", sep=""), height = 4, width = 6)
        gseaplot2(L1.YRI.HLA, geneSetID = c("L1M subfamilies", "L1P subfamilies", "L1PA subfamilies"), color = c('coral', 'red2', 'red4'), base_size = 9, pvalue_table = TRUE)
    dev.off()
    
    
    
    



    
# Save results
save(GOBP.IL16, GOBP.HLA, GOBP.HSD17B12, GOBP.SNV4, GOBP.SNV5, GOBP.YRI.HSD17B12, GOBP.YRI.HLA,
     Reactome.IL16, Reactome.HLA, Reactome.HSD17B12, Reactome.SNV4, Reactome.SNV5, Reactome.YRI.HSD17B12, Reactome.YRI.HLA,
     Hallmark.IL16, Hallmark.HLA, Hallmark.HSD17B12, Hallmark.SNV4, Hallmark.SNV5, Hallmark.YRI.HSD17B12, Hallmark.YRI.HLA,
     TFT.IL16, TFT.HLA, TFT.HSD17B12, TFT.SNV4, TFT.SNV5, TFT.YRI.HSD17B12, TFT.YRI.HLA,
     miRDB.IL16, miRDB.HLA, miRDB.HSD17B12, miRDB.SNV4, miRDB.SNV5, miRDB.YRI.HSD17B12, miRDB.YRI.HLA,
     Class.EUR.IL16, Class.EUR.HLA, Class.EUR.HSD17B12, Class.EUR.SNV4, Class.EUR.SNV5, Class.YRI.HSD17B12, Class.YRI.HLA,
     Family.EUR.IL16, Family.EUR.HLA, Family.EUR.HSD17B12, Family.EUR.SNV4, Family.EUR.SNV5, Family.YRI.HSD17B12, Family.YRI.HLA,
     L1.EUR.IL16, L1.EUR.HLA, L1.EUR.HSD17B12, L1.EUR.SNV4, L1.EUR.SNV5, L1.YRI.HSD17B12, L1.YRI.HLA,
     file = paste(my.output.dir, "All_GSEA_Results_FDR5.R", sep = ''))







# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_GSEA_Against_Top_eQTLs.txt", sep =""))
sessionInfo()
sink()      
    



    
# Clean the environment
rm(list=ls())
    
    
