# make directories
mkdir Structural_Variant_Counts
mkdir Structural_Variant_Counts/LINE1_Insertions
mkdir Structural_Variant_Counts/LINE1_deletions
mkdir Structural_Variant_Counts/ALU_Insertions
mkdir Structural_Variant_Counts/ALU_deletions


################### Extract L1/Alu Insertions and Deletions

# Go to the working directory
cd Structural_Variant_Counts

# Download Phase 3 structural variant calls and the corresponding .tbi index file
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz

# Extract autosomal TE insertions/deletions, using VCFtools
vcftools --gzvcf ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz --out AUTOSOMES.wgs.mergedSV.v8.20130502.svs.genotypes --max-missing 1 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --recode --recode-INFO-all

# Compress the output file
bgzip -@ 10 AUTOSOMES.wgs.mergedSV.v8.20130502.svs.genotypes.recode.vcf

# Extract L1 insertions
bcftools view -i 'SVTYPE ="LINE1"' AUTOSOMES.wgs.mergedSV.v8.20130502.svs.genotypes.recode.vcf.gz --threads 10 -Oz -o LINE1_Insertions/ALL_LINE1_insertions.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz

# Extract L1 deletions
bcftools view -i 'SVTYPE ="DEL_LINE1"' AUTOSOMES.wgs.mergedSV.v8.20130502.svs.genotypes.recode.vcf.gz --threads 10 -Oz -o LINE1_deletions/ALL_LINE1_deletions.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz

# Extract Alu insertions
bcftools view -i 'SVTYPE ="ALU"' AUTOSOMES.wgs.mergedSV.v8.20130502.svs.genotypes.recode.vcf.gz --threads 10 -Oz -o ALU_Insertions/ALL_Alu_insertions.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz

# Extract Alu deletions
bcftools view -i 'SVTYPE ="DEL_ALU"' AUTOSOMES.wgs.mergedSV.v8.20130502.svs.genotypes.recode.vcf.gz --threads 10 -Oz -o ALU_deletions/ALL_Alu_deletions.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz

################### Generate 0/1/2 matrices for TE insertions/deletions

# Carry out all of the commands below in the directory with the specified files

# Generate plink .bed files
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf ALL_LINE1_insertions.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz --make-bed --out ALL_LINE1_insertions.recode.plink
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf ALL_LINE1_deletions.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz --make-bed --out ALL_LINE1_deletions.recode.plink
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf ALL_Alu_insertions.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz --make-bed --out ALL_ALU_insertions.recode.plink
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf ALL_Alu_deletions.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz --make-bed --out ALL_ALU_deletions.recode.plink

# Generate .raw files, which will serve as the base for the 0/1/2 matrices
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile ALL_LINE1_insertions.recode.plink --recodeA --out ALL_LINE1_insertions.recode.plink.012
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile ALL_LINE1_deletions.recode.plink --recodeA --out ALL_LINE1_deletions.recode.plink.012
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile ALL_ALU_insertions.recode.plink --recodeA --out ALL_ALU_insertions.recode.plink.012
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile ALL_ALU_deletions.recode.plink --recodeA --out ALL_ALU_deletions.recode.plink.012

# Remove columns 2-6 from the plink.raw files
gcut --complement -d ' ' -f 2-6 ALL_LINE1_insertions.recode.plink.012.raw > ALL_LINE1_insertions.recode.plink.012
gcut --complement -d ' ' -f 2-6 ALL_LINE1_deletions.recode.plink.012.raw > ALL_LINE1_deletions.recode.plink.012
gcut --complement -d ' ' -f 2-6 ALL_ALU_insertions.recode.plink.012.raw > ALL_ALU_insertions.recode.plink.012
gcut --complement -d ' ' -f 2-6 ALL_ALU_deletions.recode.plink.012.raw > ALL_ALU_deletions.recode.plink.012

# Finally, transpose the 0/1/2 matrices
datamash -W transpose < ALL_LINE1_insertions.recode.plink.012 > ALL_LINE1_insertions_transposed.recode.plink.012
datamash -W transpose < ALL_LINE1_deletions.recode.plink.012 > ALL_LINE1_deletions_transposed.recode.plink.012
datamash -W transpose < ALL_ALU_insertions.recode.plink.012 > ALL_ALU_insertions_transposed.recode.plink.012
datamash -W transpose < ALL_ALU_deletions.recode.plink.012 > ALL_ALU_deletions_transposed.recode.plink.012


################### Run the "TE_Copy_Number_Quantification_Supplement.R" script to obtain the insertion and deletion counts
