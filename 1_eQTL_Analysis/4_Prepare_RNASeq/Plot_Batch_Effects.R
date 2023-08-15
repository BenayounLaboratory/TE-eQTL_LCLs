# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Plot_Batch_Effects_Functions.R")
    
    
    
    

   
 
# GENERAL PARAMETERS FOR PLOTS

# Define output directory
dir.output <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Quality_Control/'

# Load sample meta-data
GEUV_sample_info <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/0_All_Sample_Metadata/E-GEUV-1-unique_and_in_Phase3_snps_TE_SVs.txt", header = T, stringsAsFactors = F, sep = '\t')
rownames(GEUV_sample_info) <- GEUV_sample_info$Source.Name

# Define df for batch colors
batch_colors <- GEUV_sample_info[, c("Characteristics.ancestry.category.", "Performer", "Characteristics.sex.")]
colnames(batch_colors) <- c("ancestry.col", "lab.col", "sex.col")

# Define df for batch shapes
batch_pch <- GEUV_sample_info[, c("Characteristics.ancestry.category.", "Performer", "Characteristics.sex.")]
colnames(batch_pch) <- c("ancestry.pch", "lab.pch", "sex.pch")

# Define plot colors

    # ancestry 
    batch_colors$ancestry.col <- gsub("Tuscan", "darkorange4", batch_colors$ancestry.col)
    batch_colors$ancestry.col <- gsub("Utah", "orange4", batch_colors$ancestry.col)
    batch_colors$ancestry.col <- gsub("Finnish", "darkorange", batch_colors$ancestry.col)
    batch_colors$ancestry.col <- gsub("British", "orange", batch_colors$ancestry.col)
    batch_colors$ancestry.col <- gsub("Yoruba", "blue", batch_colors$ancestry.col)
    
    # lab 
    batch_colors$lab.col <- gsub("UNIGE", "#CC6677", batch_colors$lab.col)
    batch_colors$lab.col <- gsub("CNAG_CRG", "#117733", batch_colors$lab.col)
    batch_colors$lab.col <- gsub("MPIMG", "#88CCEE", batch_colors$lab.col)
    batch_colors$lab.col <- gsub("LUMC", "#DDCC77", batch_colors$lab.col)
    batch_colors$lab.col <- gsub("ICMB", "#332288", batch_colors$lab.col)
    batch_colors$lab.col <- gsub("HMGU", "#AA4499", batch_colors$lab.col)
    batch_colors$lab.col <- gsub("UU", "#882255", batch_colors$lab.col)
    
    # sex 
    batch_colors$sex.col <- gsub("female", "pink", batch_colors$sex.col)
    batch_colors$sex.col <- gsub("male", "dodgerblue", batch_colors$sex.col)

# Define plot shapes
    
    # ancestry
    batch_pch$ancestry.pch <- gsub("Tuscan", 18, batch_pch$ancestry.pch)
    batch_pch$ancestry.pch <- gsub("Utah", 17, batch_pch$ancestry.pch)
    batch_pch$ancestry.pch <- gsub("Finnish", 15, batch_pch$ancestry.pch)
    batch_pch$ancestry.pch <- gsub("British", 16, batch_pch$ancestry.pch)
    batch_pch$ancestry.pch <- gsub("Yoruba", 16, batch_pch$ancestry.pch)

    # lab
    batch_pch$lab.pch <- 16
    
    # sex
    batch_pch$sex.pch <- 16
    
# update plot shape values to numerics
batch_pch$ancestry.pch <- as.numeric(batch_pch$ancestry.pch)
batch_pch$lab.pch <- as.numeric(batch_pch$lab.pch)
batch_pch$sex.pch <- as.numeric(batch_pch$sex.pch)
    
    
    
    

   
 
# PLOT MAPPING RATES VS BATCHES 
      
# Parse the mapping info from each sample file
EUR_mappability <- extract_mapping_percents(dir.mapping_files = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/3_Read_Counting/Counts_EUR/Mapping_Logs_EUR/', 
                                            file.GEUV_metadata = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/0_All_Sample_Metadata/E-GEUV-1-unique_and_in_Phase3_snps_TE_SVs.txt")

YRI_mappability <- extract_mapping_percents(dir.mapping_files = '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/3_Read_Counting/Counts_YRI/Mapping_Logs_YRI/', 
                                            file.GEUV_metadata = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/0_All_Sample_Metadata/E-GEUV-1-unique_and_in_Phase3_snps_TE_SVs.txt")

# Combine mapping rate info for EUR and YRI samples
all_sample_mappability <- rbind(EUR_mappability, YRI_mappability)

# Make sure samples names are in the same order as specified by the colors df
all_sample_mappability <- all_sample_mappability[rownames(batch_colors), ]

# To improve visability, subset samples passing arbitrary mapping rate threshold 
all_sample_mappability_thresholded <- all_sample_mappability[all_sample_mappability$Unique_fraction >= 0.85, ]
  
  # Define a second color and pch df for samples passing mapping threshold
  thresholded.colors <- batch_colors[rownames(all_sample_mappability_thresholded), ]
  thresholded.pch <- batch_pch[rownames(all_sample_mappability_thresholded), ]
      
    
# Generate plots for all the data
pdf(paste(dir.output, "Plot_Mapping_Rate_vs_Batches", ".pdf", sep=""), height = 10, width = 10)
par(mfrow=c(3,3))

    # ancestry
    plot(x = all_sample_mappability$Unique_fraction, y = all_sample_mappability$Multimapping_fraction, xlab = 'Fraction Unique', ylab = 'Fraction Multimapping', main = 'Ancestry', col = batch_colors$ancestry.col, pch = batch_pch$ancestry.pch, cex = 1)
    legend('topright', legend = c('TSI', 'CEU', 'FIN', 'GBR', 'YRI'), col = c('darkorange4', 'orange4', 'darkorange', 'orange', 'blue'), pch = c(18, 17, 15, 16, 16), cex = 1, lty = 0, bty = 'n')
    
    # lab
    plot(x = all_sample_mappability$Unique_fraction, y = all_sample_mappability$Multimapping_fraction, xlab = 'Fraction Unique', ylab = 'Fraction Multimapping', main = 'Laboratory', col = batch_colors$lab.col, pch = batch_pch$lab.pch, cex = 1)
    legend('topright', legend = c('UNIGE', 'CNAG_CRG', 'MPIMG', 'LUMC', 'ICMB', 'HMGU', 'UU'), col = c("#CC6677", "#117733", "#88CCEE", "#DDCC77", "#332288", "#AA4499", "#882255"), pch = 16, cex = 1, lty = 0, bty = 'n')
    
    # Sex
    plot(x = all_sample_mappability$Unique_fraction, y = all_sample_mappability$Multimapping_fraction, xlab = 'Fraction Unique', ylab = 'Fraction Multimapping', main = 'Sex', col = batch_colors$sex.col, pch = batch_pch$sex.pch, cex = 1)
    legend('topright', legend = c('Male', 'Female'), col = c("dodgerblue", "pink"), pch = 16, cex = 1, lty = 0, bty = 'n')
    
dev.off()


# Generate plots, thresholding by mapping rate
pdf(paste(dir.output, "Plot_Mapping_Rate_vs_Batches_Zoomed", ".pdf", sep=""), height = 10, width = 10)
par(mfrow=c(3,3))

    # ancestry
    plot(x = all_sample_mappability_thresholded$Unique_fraction, y = all_sample_mappability_thresholded$Multimapping_fraction, xlab = 'Fraction Unique', ylab = 'Fraction Multimapping', main = 'Ancestry', col = thresholded.colors$ancestry.col, pch = thresholded.pch$ancestry.pch, cex = 1)
    legend('topright', legend = c('TSI', 'CEU', 'FIN', 'GBR', 'YRI'), col = c('darkorange4', 'orange4', 'darkorange', 'orange', 'blue'), pch = c(18, 17, 15, 16, 16), cex = 1, lty = 0, bty = 'n')
    
    # lab
    plot(x = all_sample_mappability_thresholded$Unique_fraction, y = all_sample_mappability_thresholded$Multimapping_fraction, xlab = 'Fraction Unique', ylab = 'Fraction Multimapping', main = 'Laboratory', col = thresholded.colors$lab.col, pch = thresholded.pch$lab.pch, cex = 1)
    legend('topright', legend = c('UNIGE', 'CNAG_CRG', 'MPIMG', 'LUMC', 'ICMB', 'HMGU', 'UU'), col = c("#CC6677", "#117733", "#88CCEE", "#DDCC77", "#332288", "#AA4499", "#882255"), pch = 16, cex = 1, lty = 0, bty = 'n')
    
    # Sex
    plot(x = all_sample_mappability_thresholded$Unique_fraction, y = all_sample_mappability_thresholded$Multimapping_fraction, xlab = 'Fraction Unique', ylab = 'Fraction Multimapping', main = 'Sex', col = thresholded.colors$sex.col, pch = thresholded.pch$sex.pch, cex = 1)
    legend('topright', legend = c('Male', 'Female'), col = c("dodgerblue", "pink"), pch = 16, cex = 1, lty = 0, bty = 'n')

dev.off()
    
    
    
    

   
 
# MDS ANALYSIS
    
# Load VST, unfiltered data
VST.YRI <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_eQTL/All_counts_YRI_86_filtered_VST.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
VST.EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/4_Prepare_RNASeq/Processed_counts_eQTL/All_counts_EUR_358_filtered_VST.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Subset gene expression
Expr.genes.YRI <- VST.YRI[grepl('ENSG', rownames(VST.YRI)), ]
Expr.genes.EUR <- VST.EUR[grepl('ENSG', rownames(VST.EUR)), ]

# Subset TE expression
Expr.TE.YRI <- VST.YRI[!grepl('ENSG', rownames(VST.YRI)) & !grepl( 'EBV', rownames(VST.YRI)), ]
Expr.TE.EUR <- VST.EUR[!grepl('ENSG', rownames(VST.EUR)) & !grepl( 'EBV', rownames(VST.EUR)), ]

# Run MDS on genes
mds.genes.YRI <- cmdscale(1-cor(Expr.genes.YRI, method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
mds.genes.EUR <- cmdscale(1-cor(Expr.genes.EUR, method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)

# Run MDS on TE data
mds.TEs.YRI <- cmdscale(1-cor(Expr.TE.YRI, method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
mds.TEs.EUR <- cmdscale(1-cor(Expr.TE.EUR, method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)

# Define MDS plot colors

    # ancestry 
    mds.col.YRI.anc <- batch_colors[rownames(mds.genes.YRI), "ancestry.col"]
    mds.col.EUR.anc <- batch_colors[rownames(mds.genes.EUR), "ancestry.col"]

    # lab 
    mds.col.YRI.lab <- batch_colors[rownames(mds.genes.YRI), "lab.col"]
    mds.col.EUR.lab <- batch_colors[rownames(mds.genes.EUR), "lab.col"]
    
    # sex 
    mds.col.YRI.sex <- batch_colors[rownames(mds.genes.YRI), "sex.col"]
    mds.col.EUR.sex <- batch_colors[rownames(mds.genes.EUR), "sex.col"]
    
# Define MDS plot point shapes

    # ancestry 
    mds.pch.YRI.anc <- batch_pch[rownames(mds.genes.YRI), "ancestry.pch"]
    mds.pch.EUR.anc <- batch_pch[rownames(mds.genes.EUR), "ancestry.pch"]

    # lab 
    mds.pch.YRI.lab <- batch_pch[rownames(mds.genes.YRI), "lab.pch"]
    mds.pch.EUR.lab <- batch_pch[rownames(mds.genes.EUR), "lab.pch"]
    
    # sex 
    mds.pch.YRI.sex <- batch_pch[rownames(mds.genes.YRI), "sex.pch"]
    mds.pch.EUR.sex <- batch_pch[rownames(mds.genes.EUR), "sex.pch"]
    
    
# Generate GENE mds plot
pdf(paste(dir.output, "Plot_MDS_genes_vs_batches", ".pdf", sep=""), width = 10, height = 10)
par(mfrow=c(3,3))

      # EUR ancestry
      plot(x = mds.genes.EUR[, 1], 
           y = mds.genes.EUR[, 2],
           pch = mds.pch.EUR.anc,
           col = mds.col.EUR.anc,
           xlab = "MDS dimension 1", 
           ylab = "MDS dimension 2",
           main = "EUR Gene Expression MDS", 
           cex = 1)
      
      legend('bottomleft', legend = c('TSI', 'CEU', 'FIN', 'GBR'), col = c('darkorange4', 'orange4', 'darkorange', 'orange'), pch = c(18, 17, 15, 16), cex = 1, lty = 0, bty = 'n')
      
      # EUR lab
      plot(x = mds.genes.EUR[, 1], 
           y = mds.genes.EUR[, 2],
           pch = mds.pch.EUR.lab,
           col = mds.col.EUR.lab,
           xlab = "MDS dimension 1", 
           ylab = "MDS dimension 2",
           main = "EUR Gene Expression MDS", 
           cex = 1)
      
      legend('bottomleft', legend = c('UNIGE', 'CNAG_CRG', 'MPIMG', 'LUMC', 'ICMB', 'HMGU', 'UU'), col = c("#CC6677", "#117733", "#88CCEE", "#DDCC77", "#332288", "#AA4499", "#882255"), pch = 16, cex = 1, lty = 0, bty = 'n')
      
      # EUR sex
      plot(x = mds.genes.EUR[, 1], 
           y = mds.genes.EUR[, 2],
           pch = mds.pch.EUR.sex,
           col = mds.col.EUR.sex,
           xlab = "MDS dimension 1", 
           ylab = "MDS dimension 2",
           main = "EUR Gene Expression MDS", 
           cex = 1)
      
      legend('bottomleft', legend = c('Male', 'Female'), col = c("blue", "pink"), pch = 16, cex = 1, lty = 0, bty = 'n')
      
      # YRI ancestry
      plot(x = mds.genes.YRI[, 1], 
           y = mds.genes.YRI[, 2],
           pch = mds.pch.YRI.anc,
           col = mds.col.YRI.anc,
           xlab = "MDS dimension 1", 
           ylab = "MDS dimension 2",
           main = "YRI Gene Expression MDS", 
           cex = 1)
      
      legend('bottomleft', legend = c('YRI'), col = c('blue'), pch = c(16), cex = 1, lty = 0, bty = 'n')
      
      # YRI lab
      plot(x = mds.genes.YRI[, 1], 
           y = mds.genes.YRI[, 2],
           pch = mds.pch.YRI.lab,
           col = mds.col.YRI.lab,
           xlab = "MDS dimension 1", 
           ylab = "MDS dimension 2",
           main = "YRI Gene Expression MDS", 
           cex = 1)
  
      legend('bottomleft', legend = c('UNIGE', 'CNAG_CRG', 'MPIMG', 'LUMC', 'ICMB', 'HMGU', 'UU'), col = c("#CC6677", "#117733", "#88CCEE", "#DDCC77", "#332288", "#AA4499", "#882255"), pch = 16, cex = 1, lty = 0, bty = 'n')
    
      # YRI sex
      plot(x = mds.genes.YRI[, 1], 
           y = mds.genes.YRI[, 2],
           pch = mds.pch.YRI.sex,
           col = mds.col.YRI.sex,
           xlab = "MDS dimension 1", 
           ylab = "MDS dimension 2",
           main = "YRI Gene Expression MDS", 
           cex = 1)
      
      legend('bottomleft', legend = c('Male', 'Female'), col = c("blue", "pink"), pch = 16, cex = 1, lty = 0, bty = 'n')
      
dev.off()
    
    
    


# Generate TE mds plot
pdf(paste(dir.output, "Plot_MDS_TEs_vs_batches", ".pdf", sep=""), width = 10, height = 10)
par(mfrow=c(3,3))

      # EUR ancestry
      plot(x = mds.TEs.EUR[, 1], 
           y = mds.TEs.EUR[, 2],
           pch = mds.pch.EUR.anc,
           col = mds.col.EUR.anc,
           xlab = "MDS dimension 1", 
           ylab = "MDS dimension 2",
           main = "EUR TE Expression MDS", 
           cex = 1)
      
      legend('bottomright', legend = c('TSI', 'CEU', 'FIN', 'GBR'), col = c('darkorange4', 'orange4', 'darkorange', 'orange'), pch = c(18, 17, 15, 16), cex = 1, lty = 0, bty = 'n')
      
      # EUR lab
      plot(x = mds.TEs.EUR[, 1], 
           y = mds.TEs.EUR[, 2],
           pch = mds.pch.EUR.lab,
           col = mds.col.EUR.lab,
           xlab = "MDS dimension 1", 
           ylab = "MDS dimension 2",
           main = "EUR TE Expression MDS", 
           cex = 1)
      
      legend('bottomright', legend = c('UNIGE', 'CNAG_CRG', 'MPIMG', 'LUMC', 'ICMB', 'HMGU', 'UU'), col = c("#CC6677", "#117733", "#88CCEE", "#DDCC77", "#332288", "#AA4499", "#882255"), pch = 16, cex = 1, lty = 0, bty = 'n')
      
      # EUR sex
      plot(x = mds.TEs.EUR[, 1], 
           y = mds.TEs.EUR[, 2],
           pch = mds.pch.EUR.sex,
           col = mds.col.EUR.sex,
           xlab = "MDS dimension 1", 
           ylab = "MDS dimension 2",
           main = "EUR TE Expression MDS", 
           cex = 1)
      
      legend('bottomright', legend = c('Male', 'Female'), col = c("blue", "pink"), pch = 16, cex = 1, lty = 0, bty = 'n')
      
      # YRI ancestry
      plot(x = mds.TEs.YRI[, 1], 
           y = mds.TEs.YRI[, 2],
           pch = mds.pch.YRI.anc,
           col = mds.col.YRI.anc,
           xlab = "MDS dimension 1", 
           ylab = "MDS dimension 2",
           main = "YRI TE Expression MDS", 
           cex = 1)
      
      legend('bottomright', legend = c('YRI'), col = c('blue'), pch = c(16), cex = 1, lty = 0, bty = 'n')
      
      # YRI lab
      plot(x = mds.TEs.YRI[, 1], 
           y = mds.TEs.YRI[, 2],
           pch = mds.pch.YRI.lab,
           col = mds.col.YRI.lab,
           xlab = "MDS dimension 1", 
           ylab = "MDS dimension 2",
           main = "YRI TE Expression MDS", 
           cex = 1)
  
      legend('bottomright', legend = c('UNIGE', 'CNAG_CRG', 'MPIMG', 'LUMC', 'ICMB', 'HMGU', 'UU'), col = c("#CC6677", "#117733", "#88CCEE", "#DDCC77", "#332288", "#AA4499", "#882255"), pch = 16, cex = 1, lty = 0, bty = 'n')
    
      # YRI sex
      plot(x = mds.TEs.YRI[, 1], 
           y = mds.TEs.YRI[, 2],
           pch = mds.pch.YRI.sex,
           col = mds.col.YRI.sex,
           xlab = "MDS dimension 1", 
           ylab = "MDS dimension 2",
           main = "YRI TE Expression MDS", 
           cex = 1)
      
      legend('bottomright', legend = c('Male', 'Female'), col = c("blue", "pink"), pch = 16, cex = 1, lty = 0, bty = 'n')
      
dev.off()


    
    
    
    
    
    
# Remove unneeded variables
rm(list=ls())
      

