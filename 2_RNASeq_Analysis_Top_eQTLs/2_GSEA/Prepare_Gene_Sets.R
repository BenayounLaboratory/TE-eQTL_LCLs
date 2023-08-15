# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(msigdbr) # Make sure this is version 7.5.1; this contains MSigDB genesets with genes in ENSEMBL format
library(stringr) # for capitalization changes
library(splitstackshape)

# Define output directory
dir.output <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/'             
            
            
          


  
  
# RETRIEVE CURATED GENE SETS  
  

# Retrieve MSigDBR HALLMARK gene sets
Hallmark_Geneset <- msigdbr(species = "Homo sapiens", category = 'H') %>% 
dplyr::select(gs_name, human_ensembl_gene)

  
# Retrieve MSigDBR REACTOME gene sets
Reactome_Geneset <- msigdbr(species = "Homo sapiens", category = 'C2', subcategory = 'CP:REACTOME') %>% 
dplyr::select(gs_name, human_ensembl_gene)
  

# Retrieve MSigDBR GO BP gene sets
GO_BP_Geneset <- msigdbr(species = "Homo sapiens", category = 'C5', subcategory = 'GO:BP') %>% 
dplyr::select(gs_name, human_ensembl_gene)


# Retrieve MSigDBR miRDB gene set
miRDB_Geneset <- msigdbr(species = "Homo sapiens", category = 'C3', subcategory = 'MIR:MIRDB') %>% 
dplyr::select(gs_name, human_ensembl_gene)
      
      
# Retrieve MSigDBR TFT GTRD gene set
GTRD_Geneset <- msigdbr(species = "Homo sapiens", category = 'C3', subcategory = 'TFT:GTRD') %>% 
dplyr::select(gs_name, human_ensembl_gene)
         






# CLEANUP GENE SET NAMES 

# remove redundant parts of labels
Hallmark_Geneset$gs_name <- gsub("HALLMARK_", "", Hallmark_Geneset$gs_name)
Reactome_Geneset$gs_name <- gsub("REACTOME_", "", Reactome_Geneset$gs_name)
GO_BP_Geneset$gs_name <- gsub("GOBP_", "", GO_BP_Geneset$gs_name)

# write function to lowercase all letters, then uppercase the first
sentence_casing <- function(input_names){
  
  # lowercase all letters
  input_lowercased <- tolower(input_names)
  
  # uppercase the first letter
  input_sentence_case <- str_to_sentence(input_lowercased)
  
  # output the final string
  return(input_sentence_case)
  
}

# change cases of gene set names, so only the first letter is capitalized
Hallmark_Geneset$gs_name <- sentence_casing(input_names = Hallmark_Geneset$gs_name)
Reactome_Geneset$gs_name <- sentence_casing(input_names = Reactome_Geneset$gs_name)
GO_BP_Geneset$gs_name <- sentence_casing(input_names = GO_BP_Geneset$gs_name)
GTRD_Geneset$gs_name <- sentence_casing(input_names = GTRD_Geneset$gs_name)

# Update underscores to either blank spaces or dashes
Hallmark_Geneset$gs_name <- gsub("_", " ", Hallmark_Geneset$gs_name)
Reactome_Geneset$gs_name <- gsub("_", " ", Reactome_Geneset$gs_name)
GO_BP_Geneset$gs_name <- gsub("_", " ", GO_BP_Geneset$gs_name)
miRDB_Geneset$gs_name <- gsub("_", "-", miRDB_Geneset$gs_name)
GTRD_Geneset$gs_name <- gsub("_", " ", GTRD_Geneset$gs_name)





  


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
  
# Make gene set (ready for GSEA)
Repeat_class.gs <- data.frame(gs_name = my.TE.classifications$Class, gene = my.TE.classifications$Full_label)
Repeat_family.gs <- data.frame(gs_name = my.TE.classifications$Family, gene = my.TE.classifications$Full_label)
           
            
          






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
  
# Make gene set (ready for GSEA)
L1_by_age.gs <- data.frame(gs_name = my.L1.classifications$Age, gene = my.L1.classifications$Full_label)

            
          



  
 
# Save gene sets
save(Hallmark_Geneset,
     Reactome_Geneset,
     GO_BP_Geneset,
     miRDB_Geneset,
     GTRD_Geneset,
     Repeat_class.gs,
     Repeat_family.gs,
     L1_by_age.gs,
     file = paste(dir.output, "Gene_Set_Collections_for_GSEA.R", sep = ''))
          




# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/2_RNASeq_Analysis_Top_eQTLs/2_GSEA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_Prepare_Gene_Sets.txt", sep =""))
sessionInfo()
sink()   

      
      
# Clean the environment
rm(list=ls())      
            


