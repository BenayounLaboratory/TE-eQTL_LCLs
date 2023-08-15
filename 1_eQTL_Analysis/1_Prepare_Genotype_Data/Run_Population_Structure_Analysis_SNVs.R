# Set strings as factors
options(stringsAsFactors = F)





# SNV POPULATION STRUCTURE USING PLINK 



# Define output directory
dir.output <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/Population_Structure_Analysis/'


# Load GEUVADIS metadata 
GEUV_sample_info <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/0_All_Sample_Metadata/E-GEUV-1-unique_and_in_Phase3_snps_TE_SVs.txt", header = T, stringsAsFactors = F, sep = '\t')
rownames(GEUV_sample_info) <- GEUV_sample_info$Source.Name

# Extract EUR and YRI metadata
EUR.samples <- GEUV_sample_info[which(GEUV_sample_info$Characteristics.ancestry.category. != 'Yoruba'), ]
YRI.samples <- GEUV_sample_info[which(GEUV_sample_info$Characteristics.ancestry.category. == 'Yoruba'), ]


# Load PCA eigenvectors.
pca.EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_PCA/Filtered_EUR_BOTH_SEXES_pca.eigenvec", header = FALSE, row.names = 1, sep = "")
pca.YRI <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_YRI_PCA/Filtered_YRI_BOTH_SEXES_pca.eigenvec", header = FALSE, row.names = 1, sep = "")

    # Match the order of the sample names
    pca.EUR <- pca.EUR[rownames(EUR.samples), ]
    pca.YRI <- pca.YRI[rownames(YRI.samples), ]

    
# Load PCA eigenvalues.
eigenval.EUR <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_PCA/Filtered_EUR_BOTH_SEXES_pca.eigenval", header = FALSE)
eigenval.YRI <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_YRI_PCA/Filtered_YRI_BOTH_SEXES_pca.eigenval", header = FALSE)


# Clean PCA df: remove sample names column and assign PC labels
pca.EUR <- pca.EUR[, -c(1)]
pca.YRI <- pca.YRI[, -c(1)]

colnames(pca.EUR) <- paste('PC', 1:ncol(pca.EUR), sep = '')
colnames(pca.YRI) <- paste('PC', 1:ncol(pca.YRI), sep = '')


# Calculate percent variance explained (pve)
pve.EUR <- 100 * (eigenval.EUR/sum(eigenval.EUR)) 
pve.YRI <- 100 * (eigenval.YRI/sum(eigenval.YRI)) 


# Add a 'color' column to the sample tables
EUR.samples$Point_color <- 'orange'
YRI.samples$Point_color <- 'blue'


    # Update the point colors
    EUR.samples[which(EUR.samples$Characteristics.ancestry.category. == 'British'), 'Point_color'] <- 'orange'
    EUR.samples[which(EUR.samples$Characteristics.ancestry.category. == 'Utah'), 'Point_color'] <- 'orange4'
    EUR.samples[which(EUR.samples$Characteristics.ancestry.category. == 'Finnish'), 'Point_color'] <- 'darkorange'
    EUR.samples[which(EUR.samples$Characteristics.ancestry.category. == 'Tuscan'), 'Point_color'] <- 'darkorange4'


# Add a 'shape' column to the sample table
EUR.samples$Point_pch <- as.numeric('16')
YRI.samples$Point_pch <- as.numeric('16')


    # Update the point shapes
    EUR.samples[which(EUR.samples$Characteristics.ancestry.category. == 'British'), 'Point_pch'] <- 16
    EUR.samples[which(EUR.samples$Characteristics.ancestry.category. == 'Utah'), 'Point_pch'] <- 17
    EUR.samples[which(EUR.samples$Characteristics.ancestry.category. == 'Finnish'), 'Point_pch'] <- 15
    EUR.samples[which(EUR.samples$Characteristics.ancestry.category. == 'Tuscan'), 'Point_pch'] <- 18



# Define plot and legend parameters
EUR.labels <- c('GBR', 'CEU', 'FIN', 'TSI')
EUR.colors <- c('orange', 'orange4', 'darkorange',  'darkorange4')
EUR.pch <- c(16, 17, 15, 18)

YRI.labels <- c('YRI')
YRI.colors <- c('blue')
YRI.pch <- c(16)

pca.label.cex <- 1
pca.axis.cex <- 1



# Plot pve
pdf(paste(dir.output, "Plot_SNV_PCA_percent_variance", ".pdf", sep=""))
par(mfrow=c(2,1))
    barplot(t(pve.EUR), names.arg = rownames(pve.EUR), xlab = "Principal component", ylab = 'Percent variance explained', main = 'EUR Genotype PCA', las = 2)
    barplot(t(pve.YRI), names.arg = rownames(pve.YRI), xlab = "Principal component", ylab = 'Percent variance explained', main = 'YRI Genotype PCA', las = 2)
dev.off()



# Plot PC vs Ancestry
pdf(paste(dir.output, "Plot_SNV_PCs_vs_Ancestry",".pdf", sep=""), width = 7, height = 7)
par(mfrow=c(2,2)) # Added to produce images of similar shape across analyses

    # EUR: PC1 vs PC2
    plot(pca.EUR$PC1,
         pca.EUR$PC2,
         col = EUR.samples$Point_color,
         pch = EUR.samples$Point_pch,
         xlab = paste('PC1', sep=""),
         ylab = paste('PC2', sep=""),
         #xlim = c(-200, 200),
         #ylim = c(-90, 90),
         cex.lab = pca.label.cex,
         cex.axis = pca.axis.cex,
         main = "EUR Population"
         )
    
        # add legend
        legend('topright', legend = EUR.labels, col = EUR.colors, pch = EUR.pch, cex = 1, lty = 0, bty = 'n')
    
    # EUR: PC1 vs PC3
    plot(pca.EUR$PC1,
         pca.EUR$PC3,
         col = EUR.samples$Point_color,
         pch = EUR.samples$Point_pch,
         xlab = paste('PC1', sep=""),
         ylab = paste('PC3', sep=""),
         #xlim = c(-200, 200),
         #ylim = c(-90, 90),
         cex.lab = pca.label.cex,
         cex.axis = pca.axis.cex,
         main = "EUR Population"
         )
    
    # YRI: PC1 vs PC2
    plot(pca.YRI$PC1,
         pca.YRI$PC2,
         col = YRI.samples$Point_color,
         pch = YRI.samples$Point_pch,
         xlab = paste('PC1', sep=""),
         ylab = paste('PC2', sep=""),
         #xlim = c(-200, 200),
         #ylim = c(-90, 90),
         cex.lab = pca.label.cex,
         cex.axis = pca.axis.cex,
         main = "YRI Population"
         )
        
        # add a legend
        legend('topright', legend = YRI.labels, col = YRI.colors, pch = YRI.pch, cex = 1, lty = 0, bty = 'n')
    
    # YRI: PC1 vs PC3
    plot(pca.YRI$PC1,
         pca.YRI$PC3,
         col = YRI.samples$Point_color,
         pch = YRI.samples$Point_pch,
         xlab = paste('PC1', sep=""),
         ylab = paste('PC3', sep=""),
         #xlim = c(-200, 200),
         #ylim = c(-90, 90),
         cex.lab = pca.label.cex,
         cex.axis = pca.axis.cex,
         main = "YRI Population"
         )
    
    
dev.off()







  
    
  
  
# Combine EUR/YRI PCs and save
pca.EUR.YRI <- rbind(pca.EUR, pca.YRI)
write.table(pca.EUR.YRI, file = paste(dir.output, "COVARIATE_SNV_Population_Structure_PCA", ".txt", sep =""), sep = "\t" , row.names = T, col.names = NA, quote=F)
  
  
  
  

  


# Clean the environment
rm(list=ls())








