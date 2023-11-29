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
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3b_RNASeq_Analysis_Candidate_Gene_OE_TElocal/3_GSEA/Comparison_repeats_for_each_treatment/'

# Load GSEA results
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3b_RNASeq_Analysis_Candidate_Gene_OE_TElocal/3_GSEA/Results_GM12878_Candidate_OE/Repeat_Families_intergenic_distal/All_GSEA_Results_FDR5_Intergenic_Distal_Repeat_Families.R')
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3b_RNASeq_Analysis_Candidate_Gene_OE_TElocal/3_GSEA/Results_GM12878_Candidate_OE/Repeat_Families_intergenic_nearby/All_GSEA_Results_FDR5_Intergenic_Nearby_Repeat_Families.R')
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3b_RNASeq_Analysis_Candidate_Gene_OE_TElocal/3_GSEA/Results_GM12878_Candidate_OE/Repeat_Families_intronic/All_GSEA_Results_FDR5_Intronic_Repeat_Families.R')
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3b_RNASeq_Analysis_Candidate_Gene_OE_TElocal/3_GSEA/Results_GM12878_Candidate_OE/Repeat_Families_exonic/All_GSEA_Results_FDR5_Exonic_Repeat_Families.R')


# Define plot labels
my.label.one <- 'Intergenic_distal'
my.label.two <- 'Intergenic_nearby'
my.label.three <- 'Intronic'
my.label.four <- 'Exonic'


# Find overlapping gene sets   
OE.IL16 <- Overlap_four_GSEA_res(result.one = distal.OE.IL16@result,
                                   result.two = nearby.OE.IL16@result,
                                   result.three = intronic.OE.IL16@result,
                                   result.four = exonic.OE.IL16@result,
                                   label.one = my.label.one,
                                   label.two = my.label.two,
                                   label.three = my.label.three,
                                   label.four = my.label.four,
                                   comparison_directionality = c('', '', '', ''),
                                   number_to_plot = 10,
                                   dir.output = my.output.dir,
                                   gs.label = "IL16_OE_Repeats")

# Find overlapping gene sets   
OE.STARD5 <- Overlap_four_GSEA_res(result.one = distal.OE.STARD5@result,
                                   result.two = nearby.OE.STARD5@result,
                                   result.three = intronic.OE.STARD5@result,
                                   result.four = exonic.OE.STARD5@result,
                                   label.one = my.label.one,
                                   label.two = my.label.two,
                                   label.three = my.label.three,
                                   label.four = my.label.four,
                                   comparison_directionality = c('', '', '', ''),
                                   number_to_plot = 10,
                                   dir.output = my.output.dir,
                                   gs.label = "STARD5_OE_Repeats")

# Find overlapping gene sets   
OE.HSD17B12 <- Overlap_four_GSEA_res(result.one = distal.OE.HSD17B12@result,
                                     result.two = nearby.OE.HSD17B12@result,
                                     result.three = intronic.OE.HSD17B12@result,
                                     result.four = exonic.OE.HSD17B12@result,
                                     label.one = my.label.one,
                                     label.two = my.label.two,
                                     label.three = my.label.three,
                                     label.four = my.label.four,
                                     comparison_directionality = c('', '', '', ''),
                                     number_to_plot = 10,
                                     dir.output = my.output.dir,
                                     gs.label = "HSD17B12_OE_Repeats")

# Find overlapping gene sets   
OE.RNF5 <- Overlap_four_GSEA_res(result.one = distal.OE.RNF5@result,
                                     result.two = nearby.OE.RNF5@result,
                                     result.three = intronic.OE.RNF5@result,
                                     result.four = exonic.OE.RNF5@result,
                                     label.one = my.label.one,
                                     label.two = my.label.two,
                                     label.three = my.label.three,
                                     label.four = my.label.four,
                                     comparison_directionality = c('', '', '', ''),
                                     number_to_plot = 10,
                                     dir.output = my.output.dir,
                                     gs.label = "RNF5_OE_Repeats")




# MAKE PLOTS

# Define general parameters

    # Dotplot point sizes
    min.dot.size <- 4
    max.dot.size <- 12


# IL16 OE 
pdf(paste(my.output.dir, "GSEA_BubblePlot_OE_IL16_Effect_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(OE.IL16, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-5,0,5)), guide = "colorbar", limits=c(-5-0.1,5+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family GSEA vs IL16 OE") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# STARD5 OE 
pdf(paste(my.output.dir, "GSEA_BubblePlot_OE_STARD5_Effect_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(OE.STARD5, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-5,0,5)), guide = "colorbar", limits=c(-5-0.1,5+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family GSEA vs STARD5 OE") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# HSD17B12 OE 
pdf(paste(my.output.dir, "GSEA_BubblePlot_OE_HSD17B12_Effect_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(OE.HSD17B12, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-5,0,5)), guide = "colorbar", limits=c(-5-0.1,5+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family GSEA vs HSD17B12 OE") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# RNF5 OE 
pdf(paste(my.output.dir, "GSEA_BubblePlot_OE_RNF5_Effect_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(OE.RNF5, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-5,0,5)), guide = "colorbar", limits=c(-5-0.1,5+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family GSEA vs RNF5 OE") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()






    
  
    
# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3b_RNASeq_Analysis_Candidate_Gene_OE_TElocal/3_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_GSEA_Comparing_Repeats.txt", sep =""))
sessionInfo()
sink()      
    





# Clean the environment
rm(list=ls())   
