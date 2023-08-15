# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Run_GSEA_Against_Top_eQTLs_Functions.R")

# Load libraries
library(ggplot2) # for bubble plot
library(scales) # for modifying the ggplot colorbar
library(reshape2) # for FC & pval bubble plot
library(DOSE)
library(poolr) # To calculate meta pvalues





# COMPARE GSEA RESULTS (rs11635336 IL16/STARD5 vs rs9271894 HLA vs rs1061810 HSD17B12) 

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Compare_Top_5_SNVs/'

# Load GSEA results
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Results_Against_Top_eQTLs/All_GSEA_Results_FDR5.R')

# Define plot labels
my.label.one <- 'rs11635336'
my.label.two <- 'rs9271894'
my.label.three <- 'rs1061810'
my.label.four <- 'rs112581165'
my.label.five <- 'rs72691418'


# Find overlapping gene sets   
overlap.Hallmark <- Overlap_five_GSEA_res(result.one = Hallmark.IL16@result,
                                         result.two = Hallmark.HLA@result,
                                         result.three = Hallmark.HSD17B12@result,
                                         result.four = Hallmark.SNV4@result,
                                         result.five = Hallmark.SNV5@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         label.three = my.label.three,
                                         label.four = my.label.four,
                                         label.five = my.label.five,
                                         comparison_directionality = c('same', 'opposite', 'opposite', 'same'),
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "Hallmark")

overlap.Reactome <- Overlap_five_GSEA_res(result.one = Reactome.IL16@result,
                                         result.two = Reactome.HLA@result,
                                         result.three = Reactome.HSD17B12@result,
                                         result.four = Reactome.SNV4@result,
                                         result.five = Reactome.SNV5@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         label.three = my.label.three,
                                         label.four = my.label.four,
                                         label.five = my.label.five,
                                         comparison_directionality = c('same', 'opposite', 'opposite', 'same'),
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "Reactome")

overlap.GOBP <- Overlap_five_GSEA_res(result.one = GOBP.IL16@result,
                                         result.two = GOBP.HLA@result,
                                         result.three = GOBP.HSD17B12@result,
                                         result.four = GOBP.SNV4@result,
                                         result.five = GOBP.SNV5@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         label.three = my.label.three,
                                         label.four = my.label.four,
                                         label.five = my.label.five,
                                         comparison_directionality = c('same', 'opposite', 'opposite', 'same'),
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "GOBP")

overlap.TFT <- Overlap_five_GSEA_res(result.one = TFT.IL16@result,
                                       result.two = TFT.HLA@result,
                                       result.three = TFT.HSD17B12@result,
                                       result.four = TFT.SNV4@result,
                                       result.five = TFT.SNV5@result,
                                       label.one = my.label.one,
                                       label.two = my.label.two,
                                       label.three = my.label.three,
                                       label.four = my.label.four,
                                       label.five = my.label.five,
                                       comparison_directionality = c('same', 'opposite', 'opposite', 'same'),
                                       number_to_plot = 10,
                                       dir.output = my.output.dir,
                                       gs.label = "TFT")

overlap.mirdb <- Overlap_five_GSEA_res(result.one = miRDB.IL16@result,
                                       result.two = miRDB.HLA@result,
                                       result.three = miRDB.HSD17B12@result,
                                       result.four = miRDB.SNV4@result,
                                       result.five = miRDB.SNV5@result,
                                       label.one = my.label.one,
                                       label.two = my.label.two,
                                       label.three = my.label.three,
                                       label.four = my.label.four,
                                       label.five = my.label.five,
                                       comparison_directionality = c('same', 'opposite', 'opposite', 'same'),
                                       number_to_plot = 10,
                                       dir.output = my.output.dir,
                                       gs.label = "miRDB")

overlap.TE_families <- Overlap_five_GSEA_res(result.one = Family.EUR.IL16@result,
                                             result.two = Family.EUR.HLA@result,
                                             result.three = Family.EUR.HSD17B12@result,
                                             result.four = Family.EUR.SNV4@result,
                                             result.five = Family.EUR.SNV5@result,
                                             label.one = my.label.one,
                                             label.two = my.label.two,
                                             label.three = my.label.three,
                                             label.four = my.label.four,
                                             label.five = my.label.five,
                                             comparison_directionality = c('same', 'opposite', 'opposite', 'same'),
                                             number_to_plot = 10,
                                             dir.output = my.output.dir,
                                             gs.label = "TEs")


  



# MAKE PLOTS

# Define general parameters

    # Dotplot point sizes
    min.dot.size <- 4
    max.dot.size <- 12


# Hallmark
pdf(paste(my.output.dir, "GSEA_BubblePlot_SNP_Effect_Overlap_Hallmark", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.Hallmark, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Hallmark Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

  
# Reactome
pdf(paste(my.output.dir, "GSEA_BubblePlot_SNP_Effect_Overlap_Reactome", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.Reactome, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Reactome Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()
  

# GOBP  
pdf(paste(my.output.dir, "GSEA_BubblePlot_SNP_Effect_Overlap_GOBP", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.GOBP, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("GO BP Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# TFT  
pdf(paste(my.output.dir, "GSEA_BubblePlot_SNP_Effect_Overlap_GTRD_TFT", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.TFT, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("TFT Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()


# miRDB
pdf(paste(my.output.dir, "GSEA_BubblePlot_SNP_Effect_Overlap_miRDB", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.mirdb, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("miRDB Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()



# TE families  
pdf(paste(my.output.dir, "GSEA_BubblePlot_SNP_Effect_Overlap_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.TE_families, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-5,0,5)), guide = "colorbar", limits=c(-5-0.1,5+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()








    
  
    
# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_GSEA_Comparing_Top_5_SNVs.txt", sep =""))
sessionInfo()
sink()      
    





# Clean the environment
rm(list=ls())   
