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
genelist.EUR.IL16 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/1_DESeq_Genotypes/DESeq_Results/EUR_rs11635336_All_Genes.txt")
genelist.EUR.HLA <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/1_DESeq_Genotypes/DESeq_Results/EUR_rs9271894_All_Genes.txt")
genelist.EUR.HSD17B12 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/1_DESeq_Genotypes/DESeq_Results/EUR_rs1061810_All_Genes.txt")
genelist.EUR.SNV4 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/1_DESeq_Genotypes/DESeq_Results/EUR_rs112581165_All_Genes.txt")
genelist.EUR.SNV5 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/1_DESeq_Genotypes/DESeq_Results/EUR_rs72691418_All_Genes.txt")
genelist.EUR.RNF5 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/1_DESeq_Genotypes/DESeq_Results/EUR_rs9270493_All_Genes.txt")
genelist.EUR.ZSCAN <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/1_DESeq_Genotypes/DESeq_Results/EUR_rs1361387_All_Genes.txt")

genelist.YRI.HSD17B12 <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/1_DESeq_Genotypes/DESeq_Results/YRI_rs2176598_All_Genes.txt")
genelist.YRI.HLA <- make_genelist(DESeq_results_path = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/1_DESeq_Genotypes/DESeq_Results/YRI_rs9271379_All_Genes.txt")



# RUN GSEA USING INTRONIC TE GENE SETS ------------------------------------------------------------------------------------------------------------

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/2_GSEA/Results_Against_Top_eQTLs/Repeat_Families_intronic/'

# Define the gene set of interest
Repeat_family.gs <- intronic_Repeat_family.gs

# Run GSEA using Repeat Family Collection
Intronic.EUR.IL16 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.IL16, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs11635336_IL16', output.object = 'unfiltered')
Intronic.EUR.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.HLA, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs9271894_HLA', output.object = 'unfiltered')
Intronic.EUR.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.HSD17B12, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs1061810_HSD17B12', output.object = 'unfiltered')
Intronic.EUR.SNV4 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.SNV4, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs112581165_SNV4', output.object = 'unfiltered')
Intronic.EUR.SNV5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.SNV5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs72691418_SNV5', output.object = 'unfiltered')
Intronic.EUR.RNF5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.RNF5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs9270493_RNF5', output.object = 'unfiltered')
Intronic.EUR.ZSCAN <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.ZSCAN, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs1361387_ZSCAN', output.object = 'unfiltered')

Intronic.YRI.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.YRI.HSD17B12, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'YRI_rs2176598_HSD17B12', output.object = 'unfiltered')
Intronic.YRI.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.YRI.HLA, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'YRI_rs9271379_HLA', output.object = 'unfiltered')

    # Generate GSEA Plots for L1 family

    # # EUR IL16 SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs11635336_IL16", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.IL16, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR HLA SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs9271894_HLA", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.HLA, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR HSD17B12 SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs1061810_HSD17B12", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.HSD17B12, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR SNV4
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs112581165_SNV4", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.SNV4, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR SNV5
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs72691418_SNV5", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.SNV5, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR RNF5 SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs9270493_RNF5", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.RNF5, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # YRI HSD17B12 SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_YRI_rs2176598_HSD17B12", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.YRI.HSD17B12, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # YRI HLA SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_YRI_rs9271379_HLA", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.YRI.HLA, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()



# Save results
save(Intronic.EUR.IL16, Intronic.EUR.HLA, Intronic.EUR.HSD17B12, Intronic.EUR.SNV4, Intronic.EUR.SNV5, Intronic.EUR.RNF5, Intronic.EUR.ZSCAN,
     Intronic.YRI.HSD17B12, Intronic.YRI.HLA,
     file = paste(my.output.dir, "GSEA_Results_FDR5_Intronic_Repeat_Families.R", sep = ''))

  
# RUN GSEA USING INTERGENIC DISTAL TE GENE SETS ------------------------------------------------------------------------------------------------------------

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/2_GSEA/Results_Against_Top_eQTLs/Repeat_Families_intergenic_distal/'

# Define the gene set of interest
Repeat_family.gs <- intergenic_distal_Repeat_family.gs

# Run GSEA using Repeat Family Collection
Distal.EUR.IL16 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.IL16, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs11635336_IL16', output.object = 'unfiltered')
Distal.EUR.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.HLA, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs9271894_HLA', output.object = 'unfiltered')
Distal.EUR.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.HSD17B12, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs1061810_HSD17B12', output.object = 'unfiltered')
Distal.EUR.SNV4 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.SNV4, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs112581165_SNV4', output.object = 'unfiltered')
Distal.EUR.SNV5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.SNV5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs72691418_SNV5', output.object = 'unfiltered')
Distal.EUR.RNF5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.RNF5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs9270493_RNF5', output.object = 'unfiltered')
Distal.EUR.ZSCAN <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.ZSCAN, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs1361387_ZSCAN', output.object = 'unfiltered')

Distal.YRI.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.YRI.HSD17B12, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'YRI_rs2176598_HSD17B12', output.object = 'unfiltered')
Distal.YRI.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.YRI.HLA, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'YRI_rs9271379_HLA', output.object = 'unfiltered')

    # Generate GSEA Plots for L1 family

    # # EUR IL16 SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs11635336_IL16", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.IL16, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR HLA SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs9271894_HLA", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.HLA, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR HSD17B12 SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs1061810_HSD17B12", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.HSD17B12, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR SNV4
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs112581165_SNV4", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.SNV4, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR SNV5
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs72691418_SNV5", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.SNV5, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR RNF5 SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs9270493_RNF5", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.RNF5, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # YRI HSD17B12 SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_YRI_rs2176598_HSD17B12", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.YRI.HSD17B12, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # YRI HLA SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_YRI_rs9271379_HLA", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.YRI.HLA, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()



# Save results
save(Distal.EUR.IL16, Distal.EUR.HLA, Distal.EUR.HSD17B12, Distal.EUR.SNV4, Distal.EUR.SNV5, Distal.EUR.RNF5, Distal.EUR.ZSCAN,
     Distal.YRI.HSD17B12, Distal.YRI.HLA,
     file = paste(my.output.dir, "GSEA_Results_FDR5_Intergenic_Distal_Repeat_Families.R", sep = ''))



# RUN GSEA USING INTERGENIC NEARBY TE GENE SETS ------------------------------------------------------------------------------------------------------------

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/2_GSEA/Results_Against_Top_eQTLs/Repeat_Families_intergenic_nearby/'

# Define the gene set of interest
Repeat_family.gs <- intergenic_nearby_Repeat_family.gs

# Run GSEA using Repeat Family Collection
Nearby.EUR.IL16 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.IL16, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs11635336_IL16', output.object = 'unfiltered')
Nearby.EUR.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.HLA, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs9271894_HLA', output.object = 'unfiltered')
Nearby.EUR.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.HSD17B12, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs1061810_HSD17B12', output.object = 'unfiltered')
Nearby.EUR.SNV4 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.SNV4, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs112581165_SNV4', output.object = 'unfiltered')
Nearby.EUR.SNV5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.SNV5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs72691418_SNV5', output.object = 'unfiltered')
Nearby.EUR.RNF5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.RNF5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs9270493_RNF5', output.object = 'unfiltered')
Nearby.EUR.ZSCAN <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.ZSCAN, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs1361387_ZSCAN', output.object = 'unfiltered')

Nearby.YRI.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.YRI.HSD17B12, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'YRI_rs2176598_HSD17B12', output.object = 'unfiltered')
Nearby.YRI.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.YRI.HLA, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'YRI_rs9271379_HLA', output.object = 'unfiltered')

    # Generate GSEA Plots for L1 family

    # # EUR IL16 SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs11635336_IL16", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.IL16, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR HLA SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs9271894_HLA", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.HLA, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR HSD17B12 SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs1061810_HSD17B12", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.HSD17B12, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR SNV4
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs112581165_SNV4", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.SNV4, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR SNV5
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs72691418_SNV5", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.SNV5, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR RNF5 SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs9270493_RNF5", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.RNF5, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # YRI HSD17B12 SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_YRI_rs2176598_HSD17B12", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.YRI.HSD17B12, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # YRI HLA SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_YRI_rs9271379_HLA", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.YRI.HLA, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()



# Save results
save(Nearby.EUR.IL16, Nearby.EUR.HLA, Nearby.EUR.HSD17B12, Nearby.EUR.SNV4, Nearby.EUR.SNV5, Nearby.EUR.RNF5, Nearby.EUR.ZSCAN,
     Nearby.YRI.HSD17B12, Nearby.YRI.HLA,
     file = paste(my.output.dir, "GSEA_Results_FDR5_Intergenic_Nearby_Repeat_Families.R", sep = ''))



# RUN GSEA USING EXONIC TE GENE SETS ------------------------------------------------------------------------------------------------------------

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/2_GSEA/Results_Against_Top_eQTLs/Repeat_Families_exonic/'

# Define the gene set of interest
Repeat_family.gs <- exonic_Repeat_family.gs

# Run GSEA using Repeat Family Collection
Exonic.EUR.IL16 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.IL16, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs11635336_IL16', output.object = 'unfiltered')
Exonic.EUR.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.HLA, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs9271894_HLA', output.object = 'unfiltered')
Exonic.EUR.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.HSD17B12, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs1061810_HSD17B12', output.object = 'unfiltered')
Exonic.EUR.SNV4 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.SNV4, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs112581165_SNV4', output.object = 'unfiltered')
Exonic.EUR.SNV5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.SNV5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs72691418_SNV5', output.object = 'unfiltered')
Exonic.EUR.RNF5 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.RNF5, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs9270493_RNF5', output.object = 'unfiltered')
Exonic.EUR.ZSCAN <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.EUR.ZSCAN, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'EUR_rs1361387_ZSCAN', output.object = 'unfiltered')

Exonic.YRI.HSD17B12 <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.YRI.HSD17B12, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'YRI_rs2176598_HSD17B12', output.object = 'unfiltered')
Exonic.YRI.HLA <- Run_GSEA(output.dir = my.output.dir, output_subfolder = '', my.genelist = genelist.YRI.HLA, my.gs.collection = Repeat_family.gs, gs.label = 'Repeat_Families', condition_label = 'YRI_rs9271379_HLA', output.object = 'unfiltered')

    # Generate GSEA Plots for L1 family

    # # EUR IL16 SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs11635336_IL16", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.IL16, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR HLA SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs9271894_HLA", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.HLA, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR HSD17B12 SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs1061810_HSD17B12", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.HSD17B12, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR SNV4
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs112581165_SNV4", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.SNV4, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR SNV5
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs72691418_SNV5", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.SNV5, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # EUR RNF5 SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_EUR_rs9270493_RNF5", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.EUR.RNF5, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # YRI HSD17B12 SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_YRI_rs2176598_HSD17B12", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.YRI.HSD17B12, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()
    # 
    # # YRI HLA SNP
    # pdf(paste(my.output.dir, "Repeat_Families/", "GSEA_EnrichmentPlot_FDR5_Repeat_Families_YRI_rs9271379_HLA", ".pdf", sep=""), height = 4, width = 6)
    #     gseaplot2(Family.YRI.HLA, geneSetID = c("L1 subfamilies"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
    # dev.off()



# Save results
save(Exonic.EUR.IL16, Exonic.EUR.HLA, Exonic.EUR.HSD17B12, Exonic.EUR.SNV4, Exonic.EUR.SNV5, Exonic.EUR.RNF5, Exonic.EUR.ZSCAN,
     Exonic.YRI.HSD17B12, Exonic.YRI.HLA,
     file = paste(my.output.dir, "GSEA_Results_FDR5_Exonic_Repeat_Families.R", sep = ''))



# SESSION INFO ------------------------------------------------------------------------------------------------------------


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/2_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_GSEA_Against_Top_eQTLs.txt", sep =""))
sessionInfo()
sink()      
    



    
# Clean the environment
rm(list=ls())
    
    
