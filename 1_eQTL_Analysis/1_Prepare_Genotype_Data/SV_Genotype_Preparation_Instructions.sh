# make directories
mkdir Structural_Variant_Genotypes
mkdir Structural_Variant_Genotypes/Plink_Pruned_EUR
mkdir Structural_Variant_Genotypes/Plink_Pruned_YRI


################### Structural Variant Genotype PCA Analysis with Plink

# Go to the working directory
cd Structural_Variant_Genotypes

# Download Phase 3 structural variant calls and the corresponding .tbi index file
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz

# Filter EUR and YRI SV genotypes using VCFTools
vcftools --vcf ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf --out INSERTIONS_FILTERED_EUR_BOTH_SEXES_chr1-22 --min-alleles 2 --max-alleles 2 --maf 0.01 --hwe .000001 --max-missing 1 --keep EUR_358_GEUV_and_Phase3_hg38 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --recode --recode-INFO-all
vcftools --vcf ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf --out INSERTIONS_FILTERED_YRI_BOTH_SEXES_chr1-22 --min-alleles 2 --max-alleles 2 --maf 0.01 --hwe .000001 --max-missing 1 --keep YRI_86_GEUV_and_Phase3_hg38 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --recode --recode-INFO-all

# Compress vcf files
bgzip -@ 10 INSERTIONS_FILTERED_EUR_BOTH_SEXES_chr1-22.recode.vcf
bgzip -@ 10 INSERTIONS_FILTERED_YRI_BOTH_SEXES_chr1-22.recode.vcf

# Prune the filtered, EUR SV genotype data
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf INSERTIONS_FILTERED_EUR_BOTH_SEXES_chr1-22.recode.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out Plink_Pruned_EUR/Insertions_Filtered_EUR_BOTH_SEXES_pruned

# Prune the filtered, YRI SV genotype data
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf INSERTIONS_FILTERED_YRI_BOTH_SEXES_chr1-22.recode.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out Plink_Pruned_YRI/Insertions_Filtered_YRI_BOTH_SEXES_pruned

# Run PCA analysis with Plink, on pruned EUR SV data
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf INSERTIONS_FILTERED_EUR_BOTH_SEXES_chr1-22.recode.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract Plink_Pruned_EUR/Insertions_Filtered_EUR_BOTH_SEXES_pruned.prune.in --pca var-wts --out Plink_Pruned_EUR/Insertions_Filtered_EUR_BOTH_SEXES_PCA

# Run PCA analysis with Plink, on pruned YRI SV data
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf INSERTIONS_FILTERED_YRI_BOTH_SEXES_chr1-22.recode.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract Plink_Pruned_YRI/Insertions_Filtered_YRI_BOTH_SEXES_pruned.prune.in --pca var-wts --out Plink_Pruned_YRI/Insertions_Filtered_YRI_BOTH_SEXES_PCA
