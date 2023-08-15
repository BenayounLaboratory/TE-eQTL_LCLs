# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Run_GSEA_Against_Top_eQTLs_Functions.R")

# Load libraries
library(ggplot2) # for bubble plot
library(scales) # for modifying the ggplot colorbar
library(reshape2) # for FC & pval bubble plot
library(DOSE)
library(poolr) # to calculate meta pvalues

# Load GSEA results
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/3_RNASeq_Analysis_Candidate_Gene_OE/4_GSEA/Results_GM12878_Candidate_OE/All_GSEA_Results_FDR5.R')
load(file = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/4_GSEA/Results_GM12878_rIL16/All_GSEA_Results_FDR5.R')


    

# COMPARE GSEA RESULTS (GM12878 rIL16 T24 vs T48) -------------------------------------------------


# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/4_GSEA/Results_Overlap/GM12878_rIL16_T24_T48/'

# Define plot labels
my.label.one <- 'rIL16_T24'
my.label.two <- 'rIL16_T48'


# Find overlapping gene sets   
overlap.GOBP <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                     result.one = GOBP.T24@result,
                                     result.two = GOBP.T48@result,
                                     label.one = my.label.one,
                                     label.two = my.label.two,
                                     number_to_plot = 10,
                                     dir.output = my.output.dir,
                                     gs.label = "GOBP")

overlap.Hallmark <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                         result.one = Hallmark.T24@result,
                                         result.two = Hallmark.T48@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "Hallmark")

overlap.Reactome <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                         result.one = Reactome.T24@result,
                                         result.two = Reactome.T48@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "Reactome")

overlap.TFT <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                   result.one = TFT.T24@result,
                                   result.two = TFT.T48@result,
                                   label.one = my.label.one,
                                   label.two = my.label.two,
                                   number_to_plot = 10,
                                   dir.output = my.output.dir,
                                   gs.label = "TFT")

overlap.miRDB <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                       result.one = miRNA.T24@result,
                                       result.two = miRNA.T48@result,
                                       label.one = my.label.one,
                                       label.two = my.label.two,
                                       number_to_plot = 10,
                                       dir.output = my.output.dir,
                                       gs.label = "miRDB")

overlap.families <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                         result.one = Family.T24@result,
                                         result.two = Family.T48@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "TEs")


    



# MAKE PLOTS

# Define general parameters

    # Dotplot point sizes
    min.dot.size <- 4
    max.dot.size <- 12
    
  
# GOBP  
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_rIL16_Overlap_GOBP", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.GOBP, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("GO BP Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# Hallmark
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_rIL16_Overlap_Hallmark", ".pdf", sep=""), width = 9, height = 5)

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
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_rIL16_Overlap_Reactome", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.Reactome, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Reactome Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()
  
# TFT
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_rIL16_Overlap_TFT", ".pdf", sep=""), width = 9, height = 5)

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
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_rIL16_Overlap_miRDB", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.miRDB, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("miRDB Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# TE families
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_rIL16_Overlap_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.families, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()







  

    





# COMPARE GSEA RESULTS (GM12878 IL16 OE vs rIL16 T24) -------------------------------------------------


# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/4_GSEA/Results_Overlap/GM12878_IL16_OE_rIL16_T24/'

# Define plot labels
my.label.one <- 'IL16_OE'
my.label.two <- 'rIL16_T24'

# Find overlapping gene sets   
overlap.GOBP <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                     result.one = GOBP.OE.IL16@result,
                                     result.two = GOBP.T24@result,
                                     label.one = my.label.one,
                                     label.two = my.label.two,
                                     number_to_plot = 10,
                                     dir.output = my.output.dir,
                                     gs.label = "GOBP")

overlap.Hallmark <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                         result.one = Hallmark.OE.IL16@result,
                                         result.two = Hallmark.T24@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "Hallmark")

overlap.Reactome <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                         result.one = Reactome.OE.IL16@result,
                                         result.two = Reactome.T24@result,
                                         label.one = my.label.one,
                                         label.two = my.label.two,
                                         number_to_plot = 10,
                                         dir.output = my.output.dir,
                                         gs.label = "Reactome")

# overlap.TFT <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
#                                     result.one = TFT.OE.IL16@result,
#                                     result.two = TFT.T24@result,
#                                     label.one = my.label.one,
#                                     label.two = my.label.two,
#                                     number_to_plot = 10)

# overlap.miRDB <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
#                                      result.one = miRDB.OE.IL16@result,
#                                      result.two = miRNA.T24@result,
#                                      label.one = my.label.one,
#                                      label.two = my.label.two,
#                                      number_to_plot = 10)

overlap.TE_families <- Overlap_two_GSEA_res(comparison_directionality = c('same'),
                                           result.one = Family.OE.IL16@result,
                                           result.two = Family.T24@result,
                                           label.one = my.label.one,
                                           label.two = my.label.two,
                                           number_to_plot = 10,
                                           dir.output = my.output.dir,
                                           gs.label = "TEs")


    



# MAKE PLOTS

# Define general parameters

    # Dotplot point sizes
    min.dot.size <- 4
    max.dot.size <- 12
    
    
# GOBP  
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_L1-positive_IL16_GOBP", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.GOBP, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("GO BP Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# Hallmark
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_L1-positive_IL16_Hallmark", ".pdf", sep=""), width = 9, height = 5)

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
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_L1-positive_IL16_Reactome", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.Reactome, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Reactome Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# # TFT
# pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_L1-positive_IL16_TFT", ".pdf", sep=""), width = 9, height = 5)
# 
#       ggplot(overlap.TFT, aes(x = Group, y = ID, colour = NES)) +
#              geom_point(aes(size = p.adjust)) +
#              scale_size(range = c(min.dot.size, max.dot.size)) +
#              scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
#              ggtitle("TFT Gene Set GSEA") +
#              labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
#              theme_bw() +
#              theme(axis.text.y = element_text(size = 9))
# 
# dev.off()

# # miRDB
# pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_L1-positive_IL16_miRDB", ".pdf", sep=""), width = 9, height = 5)
# 
#       ggplot(overlap.miRDB, aes(x = Group, y = ID, colour = NES)) +
#              geom_point(aes(size = p.adjust)) +
#              scale_size(range = c(min.dot.size, max.dot.size)) +
#              scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
#              ggtitle("miRDB Gene Set GSEA") +
#              labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
#              theme_bw() +
#              theme(axis.text.y = element_text(size = 9))
# 
# dev.off()


# TE families  
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_L1-positive_IL16_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.TE_families, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()


  

    




# COMPARE GSEA RESULTS (GM12878 All IL16) -------------------------------------------------


# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/4_GSEA/Results_Overlap/GM12878_IL16_OE_rIL16_T24_T48/'

# Define plot labels
my.label.one <- 'IL16_OE'
my.label.two <- 'rIL16_T24'
my.label.three <- 'rIL16_T48'

# Find overlapping gene sets   
overlap.GOBP <- Overlap_three_GSEA_res(comparison_directionality = c('same', 'same'),
                                       result.one = GOBP.OE.IL16@result,
                                       result.two = GOBP.T24@result,
                                       result.three = GOBP.T48@result,
                                       label.one = my.label.one,
                                       label.two = my.label.two,
                                       label.three = my.label.three,
                                       number_to_plot = 10,
                                       dir.output = my.output.dir,
                                       gs.label = "GOBP")

overlap.Hallmark <- Overlap_three_GSEA_res(comparison_directionality = c('same', 'same'),
                                           result.one = Hallmark.OE.IL16@result,
                                           result.two = Hallmark.T24@result,
                                           result.three = Hallmark.T48@result,
                                           label.one = my.label.one,
                                           label.two = my.label.two,
                                           label.three = my.label.three,
                                           number_to_plot = 10,
                                           dir.output = my.output.dir,
                                           gs.label = "Hallmark")

overlap.Reactome <- Overlap_three_GSEA_res(comparison_directionality = c('same', 'same'),
                                           result.one = Reactome.OE.IL16@result,
                                           result.two = Reactome.T24@result,
                                           result.three = Reactome.T48@result,
                                           label.one = my.label.one,
                                           label.two = my.label.two,
                                           label.three = my.label.three,
                                           number_to_plot = 10,
                                           dir.output = my.output.dir,
                                           gs.label = "Reactome")

# overlap.TFT <- Overlap_three_GSEA_res(comparison_directionality = c('same'),
#                                       result.one = TFT.OE.IL16@result,
#                                       result.two = TFT.T24@result,
#                                       label.one = my.label.one,
#                                       label.two = my.label.two,
#                                       number_to_plot = 10)

# overlap.miRDB <- Overlap_three_GSEA_res(comparison_directionality = c('same'),
#                                          result.one = miRDB.OE.IL16@result,
#                                          result.two = miRNA.T24@result,
#                                          label.one = my.label.one,
#                                          label.two = my.label.two,
#                                          number_to_plot = 10)

overlap.TE_families <- Overlap_three_GSEA_res(comparison_directionality = c('same', 'same'),
                                               result.one = Family.OE.IL16@result,
                                               result.two = Family.T24@result,
                                               result.three = Family.T48@result,
                                               label.one = my.label.one,
                                               label.two = my.label.two,
                                               label.three = my.label.three,
                                               number_to_plot = 10,
                                               dir.output = my.output.dir,
                                               gs.label = "TEs")


    



# MAKE PLOTS

# Define general parameters

    # Dotplot point sizes
    min.dot.size <- 4
    max.dot.size <- 12
    
    
# GOBP  
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_All_IL16_GOBP", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.GOBP, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("GO BP Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# Hallmark
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_All_IL16_Hallmark", ".pdf", sep=""), width = 9, height = 5)

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
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_All_IL16_Reactome", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.Reactome, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Reactome Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# # TFT
# pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_All_IL16_TFT", ".pdf", sep=""), width = 9, height = 5)
# 
#       ggplot(overlap.TFT, aes(x = Group, y = ID, colour = NES)) +
#              geom_point(aes(size = p.adjust)) +
#              scale_size(range = c(min.dot.size, max.dot.size)) +
#              scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
#              ggtitle("TFT Gene Set GSEA") +
#              labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
#              theme_bw() +
#              theme(axis.text.y = element_text(size = 9))
# 
# dev.off()

# # miRDB
# pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_All_IL16_miRDB", ".pdf", sep=""), width = 9, height = 5)
# 
#       ggplot(overlap.miRDB, aes(x = Group, y = ID, colour = NES)) +
#              geom_point(aes(size = p.adjust)) +
#              scale_size(range = c(min.dot.size, max.dot.size)) +
#              scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
#              ggtitle("miRDB Gene Set GSEA") +
#              labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
#              theme_bw() +
#              theme(axis.text.y = element_text(size = 9))
# 
# dev.off()


# TE families  
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_All_IL16_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.TE_families, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()


  

    





# COMPARE GSEA RESULTS (GM12878 IL16 OE vs STARD5 OE vs rIL16_T24) -------------------------------------------------


# GSEA directory for output files
my.output.dir <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/4_GSEA/Results_Overlap/GM12878_IL16_STARD5_OE_rIL16_T24/'

# Define plot labels
my.label.one <- 'IL16_OE'
my.label.two <- 'STARD5_OE'
my.label.three <- 'rIL16_T24'

# Find overlapping gene sets   
overlap.GOBP <- Overlap_three_GSEA_res(comparison_directionality = c('same', 'same'),
                                       result.one = GOBP.OE.IL16@result,
                                       result.two = GOBP.OE.STARD5@result,
                                       result.three = GOBP.T24@result,
                                       label.one = my.label.one,
                                       label.two = my.label.two,
                                       label.three = my.label.three,
                                       number_to_plot = 10,
                                       dir.output = my.output.dir,
                                       gs.label = "GOBP")

overlap.Hallmark <- Overlap_three_GSEA_res(comparison_directionality = c('same', 'same'),
                                           result.one = Hallmark.OE.IL16@result,
                                           result.two = Hallmark.OE.STARD5@result,
                                           result.three = Hallmark.T24@result,
                                           label.one = my.label.one,
                                           label.two = my.label.two,
                                           label.three = my.label.three,
                                           number_to_plot = 10,
                                           dir.output = my.output.dir,
                                           gs.label = "Hallmark")

overlap.Reactome <- Overlap_three_GSEA_res(comparison_directionality = c('same', 'same'),
                                           result.one = Reactome.OE.IL16@result,
                                           result.two = Reactome.OE.STARD5@result,
                                           result.three = Reactome.T24@result,
                                           label.one = my.label.one,
                                           label.two = my.label.two,
                                           label.three = my.label.three,
                                           number_to_plot = 10,
                                           dir.output = my.output.dir,
                                           gs.label = "Reactome")

# overlap.TFT <- Overlap_three_GSEA_res(comparison_directionality = c('same'),
#                                       result.one = TFT.OE.IL16@result,
#                                       result.two = TFT.T24@result,
#                                       label.one = my.label.one,
#                                       label.two = my.label.two,
#                                       number_to_plot = 10)

# overlap.miRDB <- Overlap_three_GSEA_res(comparison_directionality = c('same'),
#                                          result.one = miRDB.OE.IL16@result,
#                                          result.two = miRNA.T24@result,
#                                          label.one = my.label.one,
#                                          label.two = my.label.two,
#                                          number_to_plot = 10)

overlap.TE_families <- Overlap_three_GSEA_res(comparison_directionality = c('same', 'same'),
                                               result.one = Family.OE.IL16@result,
                                               result.two = Family.OE.STARD5@result,
                                               result.three = Family.T24@result,
                                               label.one = my.label.one,
                                               label.two = my.label.two,
                                               label.three = my.label.three,
                                               number_to_plot = 10,
                                               dir.output = my.output.dir,
                                               gs.label = "TEs")


    



# MAKE PLOTS

# Define general parameters

    # Dotplot point sizes
    min.dot.size <- 4
    max.dot.size <- 12
    
    
# GOBP  
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_All_L1_Conditions_GOBP", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.GOBP, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("GO BP Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# Hallmark
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_All_L1_Conditions_Hallmark", ".pdf", sep=""), width = 9, height = 5)

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
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_All_L1_Conditions_Reactome", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.Reactome, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Reactome Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()

# # TFT
# pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_All_L1_Conditions_TFT", ".pdf", sep=""), width = 9, height = 5)
# 
#       ggplot(overlap.TFT, aes(x = Group, y = ID, colour = NES)) +
#              geom_point(aes(size = p.adjust)) +
#              scale_size(range = c(min.dot.size, max.dot.size)) +
#              scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
#              ggtitle("TFT Gene Set GSEA") +
#              labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
#              theme_bw() +
#              theme(axis.text.y = element_text(size = 9))
# 
# dev.off()

# # miRDB
# pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_All_L1_Conditions_miRDB", ".pdf", sep=""), width = 9, height = 5)
# 
#       ggplot(overlap.miRDB, aes(x = Group, y = ID, colour = NES)) +
#              geom_point(aes(size = p.adjust)) +
#              scale_size(range = c(min.dot.size, max.dot.size)) +
#              scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
#              ggtitle("miRDB Gene Set GSEA") +
#              labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
#              theme_bw() +
#              theme(axis.text.y = element_text(size = 9))
# 
# dev.off()


# TE families  
pdf(paste(my.output.dir, "GSEA_BubblePlot_GM12878_All_L1_Conditions_TE_families", ".pdf", sep=""), width = 9, height = 5)

      ggplot(overlap.TE_families, aes(x = Group, y = ID, colour = NES)) +
             geom_point(aes(size = p.adjust)) +
             scale_size(range = c(min.dot.size, max.dot.size)) +
             scale_fill_gradientn(colours = c("blue","white","red"), values = rescale(c(-4,0,4)), guide = "colorbar", limits=c(-4-0.1,4+0.1), aesthetics = 'color') +
             ggtitle("Repeat Family Gene Set GSEA") +
             labs(x = "", y = "", colour = 'NES', size = '-log10(FDR)') +
             theme_bw() +
             theme(axis.text.y = element_text(size = 9))

dev.off()


  

    






# END CODE -------------------------------------------------


    

    
    
  
    
# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/4_RNASeq_Analysis_Peptide_Exposure/4_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_GSEA_Overlapping_Datasets.txt", sep =""))
sessionInfo()
sink()      
    



    
# Clean the environment
rm(list=ls())      
    
    
    

    





