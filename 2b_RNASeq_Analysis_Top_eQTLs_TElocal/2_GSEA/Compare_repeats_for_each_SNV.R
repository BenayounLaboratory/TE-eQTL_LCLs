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





# COMPARE GSEA RESULTS

# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/2_GSEA/Comparison_repeats_for_each_SNV/'

# Load GSEA results
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/2_GSEA/Results_Against_Top_eQTLs/Repeat_Families_intergenic_distal/GSEA_Results_FDR5_Intergenic_Distal_Repeat_Families.R')
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/2_GSEA/Results_Against_Top_eQTLs/Repeat_Families_intergenic_nearby/GSEA_Results_FDR5_Intergenic_Nearby_Repeat_Families.R')
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/2_GSEA/Results_Against_Top_eQTLs/Repeat_Families_intronic/GSEA_Results_FDR5_Intronic_Repeat_Families.R')
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/2_GSEA/Results_Against_Top_eQTLs/Repeat_Families_exonic/GSEA_Results_FDR5_Exonic_Repeat_Families.R')

# Define plot labels
my.label.one <- 'Intergenic_distal'
my.label.two <- 'Intergenic_nearby'
my.label.three <- 'Intronic'
my.label.four <- 'Exonic'


# Find overlapping gene sets   
overlap.IL16.SNV <- Overlap_four_GSEA_res(result.one = Distal.EUR.IL16@result,
                                           result.two = Nearby.EUR.IL16@result,
                                           result.three = Intronic.EUR.IL16@result,
                                           result.four = Exonic.EUR.IL16@result,
                                           label.one = my.label.one,
                                           label.two = my.label.two,
                                           label.three = my.label.three,
                                           label.four = my.label.four,
                                           comparison_directionality = c('', '', '', ''),
                                           number_to_plot = 10,
                                           dir.output = my.output.dir,
                                           gs.label = "IL16_SNV_Repeats")

overlap.HLA.SNV <- Overlap_four_GSEA_res(result.one = Distal.EUR.HLA@result,
                                           result.two = Nearby.EUR.HLA@result,
                                           result.three = Intronic.EUR.HLA@result,
                                           result.four = Exonic.EUR.HLA@result,
                                           label.one = my.label.one,
                                           label.two = my.label.two,
                                           label.three = my.label.three,
                                           label.four = my.label.four,
                                           comparison_directionality = c('', '', '', ''),
                                           number_to_plot = 10,
                                           dir.output = my.output.dir,
                                           gs.label = "HLA_SNV_Repeats")

overlap.HSD17B12.SNV <- Overlap_four_GSEA_res(result.one = Distal.EUR.HSD17B12@result,
                                               result.two = Nearby.EUR.HSD17B12@result,
                                               result.three = Intronic.EUR.HSD17B12@result,
                                               result.four = Exonic.EUR.HSD17B12@result,
                                               label.one = my.label.one,
                                               label.two = my.label.two,
                                               label.three = my.label.three,
                                               label.four = my.label.four,
                                               comparison_directionality = c('', '', '', ''),
                                               number_to_plot = 10,
                                               dir.output = my.output.dir,
                                               gs.label = "HSD17B12_SNV_Repeats")
  
overlap.SNV4 <- Overlap_four_GSEA_res(result.one = Distal.EUR.SNV4@result,
                                       result.two = Nearby.EUR.SNV4@result,
                                       result.three = Intronic.EUR.SNV4@result,
                                       result.four = Exonic.EUR.SNV4@result,
                                       label.one = my.label.one,
                                       label.two = my.label.two,
                                       label.three = my.label.three,
                                       label.four = my.label.four,
                                       comparison_directionality = c('', '', '', ''),
                                       number_to_plot = 10,
                                       dir.output = my.output.dir,
                                       gs.label = "rs112581165_SNV4_Repeats")

overlap.SNV5 <- Overlap_four_GSEA_res(result.one = Distal.EUR.SNV5@result,
                                       result.two = Nearby.EUR.SNV5@result,
                                       result.three = Intronic.EUR.SNV5@result,
                                       result.four = Exonic.EUR.SNV5@result,
                                       label.one = my.label.one,
                                       label.two = my.label.two,
                                       label.three = my.label.three,
                                       label.four = my.label.four,
                                       comparison_directionality = c('', '', '', ''),
                                       number_to_plot = 10,
                                       dir.output = my.output.dir,
                                       gs.label = "rs72691418_SNV5_Repeats")

overlap.ZSCAN <- Overlap_four_GSEA_res(result.one = Distal.EUR.ZSCAN@result,
                                       result.two = Nearby.EUR.ZSCAN@result,
                                       result.three = Intronic.EUR.ZSCAN@result,
                                       result.four = Exonic.EUR.ZSCAN@result,
                                       label.one = my.label.one,
                                       label.two = my.label.two,
                                       label.three = my.label.three,
                                       label.four = my.label.four,
                                       comparison_directionality = c('', '', '', ''),
                                       number_to_plot = 10,
                                       dir.output = my.output.dir,
                                       gs.label = "rs1361387_ZSCAN_Repeats")

overlap.YRI.HSD17B12.SNV <- Overlap_four_GSEA_res(result.one = Distal.YRI.HSD17B12@result,
                                                   result.two = Nearby.YRI.HSD17B12@result,
                                                   result.three = Intronic.YRI.HSD17B12@result,
                                                   result.four = Exonic.YRI.HSD17B12@result,
                                                   label.one = my.label.one,
                                                   label.two = my.label.two,
                                                   label.three = my.label.three,
                                                   label.four = my.label.four,
                                                   comparison_directionality = c('', '', '', ''),
                                                   number_to_plot = 10,
                                                   dir.output = my.output.dir,
                                                   gs.label = "YRI_HSD17B12_SNV_Repeats")

overlap.YRI.HLA.SNV <- Overlap_four_GSEA_res(result.one = Distal.YRI.HLA@result,
                                             result.two = Nearby.YRI.HLA@result,
                                             result.three = Intronic.YRI.HLA@result,
                                             result.four = Exonic.YRI.HLA@result,
                                             label.one = my.label.one,
                                             label.two = my.label.two,
                                             label.three = my.label.three,
                                             label.four = my.label.four,
                                             comparison_directionality = c('', '', '', ''),
                                             number_to_plot = 10,
                                             dir.output = my.output.dir,
                                             gs.label = "YRI_HLA_SNV_Repeats")

# MAKE PLOTS

# Define general parameters

    # Dotplot point sizes
    min.dot.size <- 4
    max.dot.size <- 12


# IL16 SNV
pdf(paste(my.output.dir, "GSEA_BubblePlot_IL16_SNV_Effect_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.IL16.SNV, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-5,0,5)), guide = "colorbar", limits=c(-5-0.1,5+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family GSEA vs IL16/STARD5 SNV") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# HLA SNV
pdf(paste(my.output.dir, "GSEA_BubblePlot_HLA_SNV_Effect_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.HLA.SNV, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-5,0,5)), guide = "colorbar", limits=c(-5-0.1,5+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family GSEA vs HLA SNV") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# HSD17B12 SNV
pdf(paste(my.output.dir, "GSEA_BubblePlot_HSD17B12_SNV_Effect_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.HSD17B12.SNV, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-5,0,5)), guide = "colorbar", limits=c(-5-0.1,5+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family GSEA vs HSD17B12 SNV") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# rs112581165 SNV4
pdf(paste(my.output.dir, "GSEA_BubblePlot_rs112581165_SNV4_Effect_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.SNV4, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-5,0,5)), guide = "colorbar", limits=c(-5-0.1,5+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family GSEA vs rs112581165 SNV4") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# rs72691418 SNV5
pdf(paste(my.output.dir, "GSEA_BubblePlot_rs72691418_SNV5_Effect_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.SNV5, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-5,0,5)), guide = "colorbar", limits=c(-5-0.1,5+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family GSEA vs rs72691418 SNV5") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# rs1361387 ZSCAN
pdf(paste(my.output.dir, "GSEA_BubblePlot_rs1361387_ZSCAN_Effect_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.ZSCAN, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-5,0,5)), guide = "colorbar", limits=c(-5-0.1,5+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family GSEA vs rs1361387 ZSCAN") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# YRI HSD17B12 SNV
pdf(paste(my.output.dir, "GSEA_BubblePlot_YRI_HSD17B12_SNV_Effect_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.YRI.HSD17B12.SNV, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-5,0,5)), guide = "colorbar", limits=c(-5-0.1,5+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family GSEA vs HSD17B12 SNV") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# YRI HLA SNV
pdf(paste(my.output.dir, "GSEA_BubblePlot_YRI_HLA_SNV_Effect_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.YRI.HLA.SNV, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-5,0,5)), guide = "colorbar", limits=c(-5-0.1,5+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family GSEA vs HLA SNV") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()


    
  
    
# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs_TElocal/2_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_GSEA_Comparing_Repeats.txt", sep =""))
sessionInfo()
sink()      
    





# Clean the environment
rm(list=ls())   
