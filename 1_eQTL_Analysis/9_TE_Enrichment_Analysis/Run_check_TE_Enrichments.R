# Set strings as factors
options(stringsAsFactors = F)

# Define output directory
dir.output <- '/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/9_TE_Enrichment_Analysis/Results/'


# PREPARE BED FILES FOR TESTING (RUN THESE COMMANDS IN TERMINAL) #################################################################################################################################################################################################################################

# Convert the Hammel lab TE GTF file to a bed file
# gtf2bed < GRCh38_GENCODE_rmsk_TE_2020_Feb_19.gtf > GRCh38_GENCODE_rmsk_TE_2020_Feb_19.bed

# Make BED file for polymorphic structural variants 
# cat ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf |tail -n +12 |awk '{FS="\t";OFS="\t";print $1,$2-1,$2,$3, etc}' > Structural_Variants_1000G.bed

    # afterwards, manually delete the header and the last empty line





# PREPARE SNP FILES FOR TESTING #################################################################################################################################################################################################################################


# Read in ALL significant trans-eQTL SNPs found by the lab
sig_snps <- read.csv(file="/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/6_Integrate_cis_trans_eQTLs/Annotated_and_filtered_eQTLs/EUR_L1_eQTLs_annotated_0_filters_2023-05-22.txt", header=TRUE, sep='\t') 
    
    # Keep only unique RsIDs
    sig_snps <- unique(sig_snps$RsID)
    
# Read in all SNVs tested in the eQTL analysis
all_snps <- read.csv(file="/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/RESOURCE_SNP_info_EUR.txt", header = TRUE, sep = '\t')

    # Update colnames
    colnames(all_snps)[which(names(all_snps) == "RsID")] <- "rs_id"
    
    # Assign RsID to rownames
    rownames(all_snps) <- all_snps$rs_id

# Read in TOP significant trans-eQTL index SNPs
top.sig_snps <- read.csv(file="/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/5_Scan_for_eQTLs/Significant_SNVs/L1_trans_eQTLs_SIGNIFICANT_EUR.clumped", header = TRUE, row.names = NULL, stringsAsFactors = F, sep = "") 

# Generate BED files for TOP sig SNPs
top.sig_snps.bed <- all_snps[unique(top.sig_snps$SNP), ]

    # Keep bed columns
    top.sig_snps.bed <- top.sig_snps.bed[, c('chr', 'pos', 'pos')]
    
    # Add name
    top.sig_snps.bed$Name <- rownames(top.sig_snps.bed)
    
    # Update colnames
    colnames(top.sig_snps.bed) <- c('chr', 'start', 'end', 'name')
    
    # Add 10,000 bp windows (5000 bp on each side)
    top.sig_snps.bed$start <- top.sig_snps.bed$start - 5000
    top.sig_snps.bed$end <- top.sig_snps.bed$end + 5000
    
        # Save as BED file, noting that chr names don't have the chr label
        write.table(top.sig_snps.bed, file = paste(dir.output, 'My_Index_SNVs_10000_window_without_chr', '.bed', sep =""), sep = "\t", row.names = F, col.names = FALSE, quote = F)
        
    # Add chr label to chromosomes
    top.sig_snps.bed$chr <- paste('chr', top.sig_snps.bed$chr, sep = '')
    
        # Save as BED file, noting that chr names have the chr label
        write.table(top.sig_snps.bed, file = paste(dir.output, 'My_Index_SNVs_10000_window_with_chr', '.bed', sep =""), sep = "\t", row.names = F, col.names = FALSE, quote = F)
    
    
# Generate BED files for RANDOM SNPs
    
    # Remove eQTL SNPs from the SNP list
    all_snps <- all_snps[!(all_snps$rs_id %in% sig_snps), ]
    
    # Set set for reproducibility of randomness
    set.seed(90280)
    
    # Find 1000 random indices
    my.random.indices <- runif(n = 1000, min = 1, max = nrow(all_snps))
    
    # Extract SNP info
    random.bed <- all_snps[my.random.indices, ]
    
    # Keep bed columns
    random.bed <- random.bed[, c('chr', 'pos', 'pos')]
    
    # Add name column (RsID)
    random.bed$Name <- rownames(random.bed)
    
    # Update colnames
    colnames(random.bed) <- c('chr', 'start', 'end', 'name')
    
    # Add 10,000 bp windows (5000 bp on each side)
    random.bed$start <- random.bed$start - 5000
    random.bed$end <- random.bed$end + 5000

        # Save as BED file, noting that chr names don't have the chr label
        write.table(random.bed, file = paste(dir.output, 'My_Random_SNVs_10000_window_without_chr', '.bed', sep =""), sep = "\t", row.names = F, col.names = FALSE, quote = F)
            
    # Add chr label to chromosomes
    random.bed$chr <- paste('chr', random.bed$chr, sep = '')
    
        # Save as BED file, noting that chr names have the chr label
        write.table(random.bed, file = paste(dir.output, 'My_Random_SNVs_10000_window_with_chr', '.bed', sep =""), sep = "\t", row.names = F, col.names = FALSE, quote = F)
    

# Remove large variables from memory
 rm(all_snps)       
        
 
 
 
        
# CHECK FOR SNP WINDOW AND TE OVERLAPS (RUN THESE COMMANDS IN TERMINAL) #################################################################################################################################################################################################################################
   
# Overlap TOP eQTLs with Hammel lab TE annotations
#bedtools intersect -a My_Index_SNVs_10000_window_with_chr.bed -b /Users/juanb/Library/CloudStorage/Dropbox/Research/Shared_Bioinformatic_Resources/STAR_genome_indices/TEtranscripts_rmsk_files/GRCh38_GENCODE_rmsk_TE_2020_Feb_19.bed -wb > Results_Top_eQTLs_overlapping_fixed_TEs.txt

# Overlap RANDOM SNVs with Hammel lab TE annotations
#bedtools intersect -a My_Random_SNVs_10000_window_with_chr.bed -b /Users/juanb/Library/CloudStorage/Dropbox/Research/Shared_Bioinformatic_Resources/STAR_genome_indices/TEtranscripts_rmsk_files/GRCh38_GENCODE_rmsk_TE_2020_Feb_19.bed -wb > Results_Random_SNVs_overlapping_fixed_TEs.txt
    
# Overlap TOP eQTLs with 1000G structural variants
#bedtools intersect -a My_Index_SNVs_10000_window_without_chr.bed -b Structural_Variants_1000G.bed -wb > Results_Top_eQTLs_overlapping_polymorphic_SVs.txt

# Overlap RANDOM SNVs with 1000G structural variants
#bedtools intersect -a My_Random_SNVs_10000_window_without_chr.bed -b Structural_Variants_1000G.bed -wb > Results_Random_SNVs_overlapping_polymorphic_SVs.txt

# NOTE: L1 structural variants that overlap either my eQTLs or random SNVs are limited. For simplicity, I'll only focus on the reference TEs below.

 
 
# PLOT DATA #################################################################################################################################################################################################################################

 
# Load SNP/TE overlaps
overlaps.eQTLs.fixedTEs <- read.csv(file="/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/9_TE_Enrichment_Analysis/Results/Results_Top_eQTLs_overlapping_fixed_TEs.txt", header = FALSE, sep='\t') 
overlaps.Random.fixedTEs <- read.csv(file="/Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/9_TE_Enrichment_Analysis/Results/Results_Random_SNVs_overlapping_fixed_TEs.txt", header = FALSE, sep='\t') 

# Subset L1 overlaps
overlaps.eQTLs.fixedTEs <- overlaps.eQTLs.fixedTEs[grepl("family_id L1", overlaps.eQTLs.fixedTEs$V14), ]
overlaps.Random.fixedTEs <- overlaps.Random.fixedTEs[grepl("family_id L1", overlaps.Random.fixedTEs$V14), ]

# Add column to count the number L1 fragments per SNV
overlaps.eQTLs.fixedTEs$TE_count <- 1
overlaps.Random.fixedTEs$TE_count <- 1

# Subset only the SNVs and counts
overlaps.eQTLs.fixedTEs.counts <- overlaps.eQTLs.fixedTEs[, c('V4', 'TE_count')]
overlaps.Random.fixedTEs.counts <- overlaps.Random.fixedTEs[, c('V4', 'TE_count')]

# Aggregate the overlaps by SNV, adding up the counts
overlaps.eQTLs.fixedTEs.counts <- aggregate(overlaps.eQTLs.fixedTEs.counts[, -c(1)], 
                                            by = list(overlaps.eQTLs.fixedTEs.counts$V4), 
                                            FUN = 'sum')

overlaps.Random.fixedTEs.counts <- aggregate(overlaps.Random.fixedTEs.counts[, -c(1)], 
                                            by = list(overlaps.Random.fixedTEs.counts$V4), 
                                            FUN = 'sum')

# Update colnames
colnames(overlaps.eQTLs.fixedTEs.counts) <- c('RsID', 'TE_count')
colnames(overlaps.Random.fixedTEs.counts) <- c('RsID', 'TE_count')

# Assign RsID to rownames
rownames(overlaps.eQTLs.fixedTEs.counts) <- overlaps.eQTLs.fixedTEs.counts$RsID
rownames(overlaps.Random.fixedTEs.counts) <- overlaps.Random.fixedTEs.counts$RsID

# Add a coloumn to the BED files to hold the final counts (since SNVs that have 0 overlap are not reprented in the current counts table)
top.sig_snps.bed$TE_count <- 0
random.bed$TE_count <- 0

# Add a coloumn to the BED files to hold whether results are from Random or eQTLs (since SNVs that have 0 overlap are not reprented in the current counts table)
top.sig_snps.bed$Group <- 'Index_SNVs'
random.bed$Group <- 'Random_SNVs'

# Transfer counts to the BED table
top.sig_snps.bed[rownames(overlaps.eQTLs.fixedTEs.counts), 'TE_count'] <- overlaps.eQTLs.fixedTEs.counts$TE_count
random.bed[rownames(overlaps.Random.fixedTEs.counts), 'TE_count'] <- overlaps.Random.fixedTEs.counts$TE_count

# Combine the counts from the random and real data
all_overlap_counts <- rbind(top.sig_snps.bed[, c('Group', 'TE_count')], random.bed[, c('Group', 'TE_count')])

# Reorder the group levels in the dataframe to match the order I want in the plot
group_order <- c('Random_SNVs', 'Index_SNVs') 
all_overlap_counts$Group <- factor(all_overlap_counts$Group , levels = group_order)

# Run statistical analyses
my.stats <- wilcox.test(all_overlap_counts[which(all_overlap_counts$Group == 'Index_SNVs'), 'TE_count'],
                        all_overlap_counts[which(all_overlap_counts$Group == 'Random_SNVs'), 'TE_count'])

# Define general plots parameters

    # Specify x tick positions
    plot.other.ticks <- c(1, 2)
    
    # Define point colors
    my.main.colors <- c('black', 'red')
    other.point.color <- rep(my.main.colors, 1000)
    
    # Define groups
    plot.other.label <- c('Random_SNVs', 'Index_SNVs')
    
# Start PDF to plot other responses
pdf(paste(dir.output, "Boxplot_TE_Overlap_Counts", ".pdf", sep=""), width = 8, height = 8)
#par(mfrow = c(2, 2)) # 2x2 plotting matrix

    # make boxplots 
    boxplot(TE_count ~ Group, 
            data = all_overlap_counts,
            ylim = c(0, 25),
            outline = F, 
            #col = boxplot.col,
            #main = paste('RT-qPCR Expression', sep = ''),
            ylab = c('Number of L1 Fragments within +/- 5 kb'),
            xlab = c(''),
            xaxt = 'n')
    
    # overlay the points
    stripchart(TE_count ~ Group, data = all_overlap_counts, vertical = TRUE, add = TRUE, method = 'jitter', pch = 20, cex = 1, col = other.point.color) 
    
    # Specify axis labels
    axis(1, at = plot.other.ticks, labels = plot.other.label)
    
    # Add pvalues
    text(1.5, 24.5, signif(my.stats$p.value, 3), cex = 1)
    
    
# End pdf
dev.off()

