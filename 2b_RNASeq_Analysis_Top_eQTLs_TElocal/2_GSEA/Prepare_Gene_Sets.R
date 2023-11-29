# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(stringr) # for capitalization changes
library(splitstackshape)

# Define output directory
dir.output <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs/2_GSEA/'             
        
  
  






  


# MAKE TRANSPOSON GENE SETS (BY CLASS AND FAMILY, USING TETRANSCRIPTS DESIGNATIONS)

# Load TETranscripts raw counts file (which has all the annotated transposons)  
my.counts <- read.csv(file = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/3_Read_Counting/Counts_EUR/fastp_ERR188021Aligned.sortedByCoord.out.bam.cntTable", header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = '\t')
  
# Define rows that are TEs
TE.rows <- which(!grepl('ENSG', rownames(my.counts)) & !grepl('EBV', rownames(my.counts)) & !grepl('pcDNA', rownames(my.counts)))

# Make dataframe to hold TE classifications
my.TE.classifications <- data.frame(row.names = rownames(my.counts[TE.rows, , drop = FALSE]),
                                    Full_label = rownames(my.counts[TE.rows, , drop = FALSE]),
                                    To_split = rownames(my.counts[TE.rows, , drop = FALSE]))
  
# Split TE labels into subfamily, family, and class. 
my.TE.classifications <- as.data.frame(splitstackshape::cSplit(my.TE.classifications, splitCols = "To_split", sep = ":", direction = "wide", fixed = TRUE, drop = TRUE, stripWhite = FALSE, makeEqual = TRUE))

# Update the column names
names(my.TE.classifications)[2:4] <- c('Subfamily', 'Family', 'Class')   

# Append "subfamilies" to the Family/Class names, since these names will be used as the gene set names
my.TE.classifications$Class <- paste(my.TE.classifications$Class, ' subfamilies', sep = '')
my.TE.classifications$Family <- paste(my.TE.classifications$Family, ' subfamilies', sep = '')

# Make general gene set 
Repeat_class.gs <- data.frame(gs_name = my.TE.classifications$Class, gene = my.TE.classifications$Full_label)
Repeat_family.gs <- data.frame(gs_name = my.TE.classifications$Family, gene = my.TE.classifications$Full_label)
           
# Make dummy gene sets for various TE loci
intergenic_distal_Repeat_family.gs <- Repeat_family.gs
intergenic_nearby_Repeat_family.gs <- Repeat_family.gs
intronic_Repeat_family.gs <- Repeat_family.gs
exonic_Repeat_family.gs <- Repeat_family.gs

    # Update the TE labels
    intergenic_distal_Repeat_family.gs$gene <- paste('intergenic_distal', intergenic_distal_Repeat_family.gs$gene, sep = ":")
    intergenic_nearby_Repeat_family.gs$gene <- paste('intergenic_near', intergenic_nearby_Repeat_family.gs$gene, sep = ":")
    intronic_Repeat_family.gs$gene <- paste('intronic', intronic_Repeat_family.gs$gene, sep = ":")
    exonic_Repeat_family.gs$gene <- paste('exonic', exonic_Repeat_family.gs$gene, sep = ":")



    
    
# MAKE L1 GENE SETS (BY L1 AGE, USING TETRANSCRIPTS DESIGNATIONS)
  
# Load TETranscripts raw counts file (which has all the annotated transposons)  
my.counts <- read.csv(file = "/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/3_Read_Counting/Counts_EUR/fastp_ERR188021Aligned.sortedByCoord.out.bam.cntTable", header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = '\t')
  
# Define rows that are L1
L1.rows <- which(grepl(':L1:', rownames(my.counts)))

# Make dataframe to hold TE classifications
my.L1.classifications <- data.frame(row.names = rownames(my.counts[L1.rows, , drop = FALSE]),
                                    Full_label = rownames(my.counts[L1.rows, , drop = FALSE]),
                                    To_split = rownames(my.counts[L1.rows, , drop = FALSE]))


# Segregate by age
L1M_subfamilies <- my.L1.classifications[which(grepl('L1M', rownames(my.L1.classifications))), ]
L1P_subfamilies <- my.L1.classifications[which(grepl('L1P', rownames(my.L1.classifications)) & !grepl('L1PA', rownames(my.L1.classifications))), ]
L1PA_subfamilies <- my.L1.classifications[which(grepl('L1PA', rownames(my.L1.classifications)) | grepl('L1H', rownames(my.L1.classifications))), ]

# Add age columns
L1M_subfamilies$Age <- 'L1M subfamilies'
L1P_subfamilies$Age <- 'L1P subfamilies'
L1PA_subfamilies$Age <- 'L1PA subfamilies'

# Rejoin L1 subfamilies
my.L1.classifications <- rbind(L1M_subfamilies, L1P_subfamilies, L1PA_subfamilies)
  
# Make a general gene set
L1_by_age.gs <- data.frame(gs_name = my.L1.classifications$Age, gene = my.L1.classifications$Full_label)

# Make dummy gene sets for various TE loci
intergenic_distal_L1_by_age.gs <- L1_by_age.gs
intergenic_nearby_L1_by_age.gs <- L1_by_age.gs
intronic_L1_by_age.gs <- L1_by_age.gs
exonic_L1_by_age.gs <- L1_by_age.gs

    # Update the TE labels
    intergenic_distal_L1_by_age.gs$gene <- paste('intergenic_distal', intergenic_distal_L1_by_age.gs$gene, sep = ":")
    intergenic_nearby_L1_by_age.gs$gene <- paste('intergenic_near', intergenic_nearby_L1_by_age.gs$gene, sep = ":")
    intronic_L1_by_age.gs$gene <- paste('intronic', intronic_L1_by_age.gs$gene, sep = ":")
    exonic_L1_by_age.gs$gene <- paste('exonic', exonic_L1_by_age.gs$gene, sep = ":")         
          



  
 
# Save gene sets
save(intergenic_distal_Repeat_family.gs,
     intergenic_nearby_Repeat_family.gs,
     intronic_Repeat_family.gs,
     exonic_Repeat_family.gs,
     intergenic_distal_L1_by_age.gs,
     intergenic_nearby_L1_by_age.gs,
     intronic_L1_by_age.gs,
     exonic_L1_by_age.gs,
     file = paste(dir.output, "Repeat_Subset_Gene_Set_Collections_for_GSEA.R", sep = ''))
          




# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2b_RNASeq_Analysis_Top_eQTLs/2_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_Prepare_Gene_Sets.txt", sep =""))
sessionInfo()
sink()   

      
      
# Clean the environment
rm(list=ls())      
            


