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
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/3_GSEA/Comparison_repeats_for_each_treatment/'

# Load GSEA results
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/3_GSEA/Results_GM12878_rIL16/Repeat_Families_intergenic_distal/All_GSEA_Results_FDR5_Intergenic_Distal_Repeat_Families.R')
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/3_GSEA/Results_GM12878_rIL16/Repeat_Families_intergenic_nearby/All_GSEA_Results_FDR5_Intergenic_Nearby_Repeat_Families.R')
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/3_GSEA/Results_GM12878_rIL16/Repeat_Families_intronic/All_GSEA_Results_FDR5_Intronic_Repeat_Families.R')
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/3_GSEA/Results_GM12878_rIL16/Repeat_Families_exonic/All_GSEA_Results_FDR5_Exonic_Repeat_Families.R')


# Define plot labels
my.label.one <- 'Intergenic_distal'
my.label.two <- 'Intergenic_nearby'
my.label.three <- 'Intronic'
my.label.four <- 'Exonic'


# Find overlapping gene sets   
T24.rIL16 <- Overlap_four_GSEA_res(result.one = distal.T24@result,
                                   result.two = nearby.T24@result,
                                   result.three = intronic.T24@result,
                                   result.four = exonic.T24@result,
                                   label.one = my.label.one,
                                   label.two = my.label.two,
                                   label.three = my.label.three,
                                   label.four = my.label.four,
                                   comparison_directionality = c('', '', '', ''),
                                   number_to_plot = 10,
                                   dir.output = my.output.dir,
                                   gs.label = "rIL16_T24_Repeats")

T48.rIL16 <- Overlap_four_GSEA_res(result.one = distal.T48@result,
                                   result.two = nearby.T48@result,
                                   result.three = intronic.T48@result,
                                   result.four = exonic.T48@result,
                                   label.one = my.label.one,
                                   label.two = my.label.two,
                                   label.three = my.label.three,
                                   label.four = my.label.four,
                                   comparison_directionality = c('', '', '', ''),
                                   number_to_plot = 10,
                                   dir.output = my.output.dir,
                                   gs.label = "rIL16_T48_Repeats")


# MAKE PLOTS

# Define general parameters

    # Dotplot point sizes
    min.dot.size <- 4
    max.dot.size <- 12


# T24 rIL16 
pdf(paste(my.output.dir, "GSEA_BubblePlot_rIL16_T24_Effect_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(T24.rIL16, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-5,0,5)), guide = "colorbar", limits=c(-5-0.1,5+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family GSEA vs rIL16 T24") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# T48 rIL16 
pdf(paste(my.output.dir, "GSEA_BubblePlot_rIL16_T48_Effect_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(T48.rIL16, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-5,0,5)), guide = "colorbar", limits=c(-5-0.1,5+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family GSEA vs rIL16 T48") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()




    
  
    
# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4b_RNASeq_Analysis_Peptide_Exposure_TElocal/3_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_GSEA_Comparing_Repeats.txt", sep =""))
sessionInfo()
sink()      
    





# Clean the environment
rm(list=ls())   
