# Set strings as factors
options(stringsAsFactors = F)

# Define the output directory
dir.output <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/5_GWAS_Analyses/2_Visualize_Collaborator_Results/Plots/'

# Load Suh Lab GWAS summary results
GWAS_summary_res <- read.csv("/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/5_GWAS_Analyses/2_Visualize_Collaborator_Results/Suh_Lab_Results/Suh_Lab_Results_Modified.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Subset data with an age-trait association
age.trait.summary.res <- GWAS_summary_res[which(GWAS_summary_res$Have.age.related.triats == c('Yes')), ]





# PIE CHART FOR FRACTION OF SNPS WITH AGE-RELATED TRAIT

# Count the number of snps
number_no_mesh <- length(which(GWAS_summary_res$Have.age.related.triats == c('No')))
number_yes_mesh <- length(which(GWAS_summary_res$Have.age.related.triats == c('Yes')))

# Collect the numbers in one vector
has_mesh_ID <- c(No = number_no_mesh, Yes = number_yes_mesh)

# Calculate percents of total
percent_no_mesh <- 100 * (number_no_mesh / length(GWAS_summary_res$Have.age.related.triats))
percent_yes_mesh <- 100 * (number_yes_mesh / length(GWAS_summary_res$Have.age.related.triats))

# Keep only 1 decimal
percent_no_mesh <- signif(percent_no_mesh, digits = 3)
percent_yes_mesh <- signif(percent_yes_mesh, digits = 3)

# Define plot labels
my_pie_labels <- c(paste('No ', percent_no_mesh, '%', sep = ''),
                paste('Yes ', percent_yes_mesh, '%', sep = ''))

# Make a pie chart for the # of SNPs with or without MeSH ID
pdf(paste(dir.output, "Pie_Chart_Percent_of_SNPs_Age-Related_Traits",".pdf", sep=""), width = 4, height = 4)

pie(has_mesh_ID, 
    col = c('#FFFFFF', '#808080'), 
    labels = my_pie_labels,
    main = c("SNVs with Age-Related Traits"))

dev.off()  





# HISTOGRAM WITH # PHEWAS MESH IDS 
pdf(paste(dir.output, "Histogram_#_of_PheWAS_MeSH_IDs",".pdf", sep=""), width = 5, height = 5)

hist(age.trait.summary.res$N.of.PheWAS.MeSH_IDs,
     main = '# of Associated Age-Related Traits per SNV',
     ylab = '# of SNVs',
     ylim = c(0, 120),
     xlab = '# of PheWAS MeSH IDs',
     xlim = c(0, 35),
     breaks = 35,
     col = '#808080',
     freq = TRUE)

dev.off()  





# HISTOGRAM WITH # OF LINKED MESH IDS
pdf(paste(dir.output, "Histogram_#_of_Linked_MeSH_IDs",".pdf", sep=""), width = 5, height = 5)

hist(age.trait.summary.res$No.of.Linked.MeSH_IDs,
     main = '# of Linked Age-Related Traits per SNV',
     ylab = '# of SNVs',
     ylim = c(0, 120),
     xlab = '# of Linked MeSH IDs',
     xlim = c(0, 5),
     breaks = 4,
     col = '#808080',
     freq = TRUE)

dev.off()  

    
  

    



    
# Clean the environment
rm(list=ls())






