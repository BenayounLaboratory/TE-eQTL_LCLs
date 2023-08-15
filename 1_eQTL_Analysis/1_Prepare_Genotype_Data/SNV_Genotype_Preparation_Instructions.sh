# make directories
mkdir SNV_Genotypes
mkdir SNV_Genotypes/Plink_EUR_Genotypes
mkdir SNV_Genotypes/Plink_EUR_PCA
mkdir SNV_Genotypes/Plink_YRI_Genotypes
mkdir SNV_Genotypes/Plink_YRI_PCA
mkdir SNV_Genotypes/dbSNP_Annotations


################### Prepare dbSNP build 155

# Go to directory where dbSNP data will be processed
cd SNV_Genotypes/dbSNP_Annotations

# Download the dbSNP VCF file (GCF_000001405.39.gz was used for this study), as well as the corresponding .tbi index file
# Link: https://ftp.ncbi.nih.gov/snp/latest_release/VCF/

# Download the corresponding genome assembly report (this will be used to map Refseq IDs in the dbSNP file to chromosome numbers)
# Link: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt

# Make a file with the chromosome refseq ids and chromosome number mappings
for k in *assembly_report.txt
  do
    out=$(echo $k | sed 's/.txt/.chrnames/')
    grep -e '^[^#]' $k | awk '{ print $7, $1 }' > $out
done

# Update the dbSNP file with chromosome numbers
bcftools annotate --rename-chrs GCF_000001405.39_GRCh38.p13_assembly_report.chrnames --threads 10 -Oz -o GRCh38.dbSNP155.vcf.gz GCF_000001405.39.gz


################### Annotate 1000Genomes SNVs with RsIDs

# Go up to the main genotype directory
cd ..

# Download GRCh38 1000Genomes Phase 3 SNV/Indel Data
# Link: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/supporting/related_samples/

# Index SNV files
for f in *.vcf.gz;do bcftools index $f --threads 10;done

# Apply dbSNP155 RsIDs
for f in *.vcf.gz;do bcftools annotate -a dbSNP_Annotations/GRCh38.dbSNP155.vcf.gz -c ID $f --threads 10 -Oz -o rsID.${f};done

# Delete or move the original SNV files elsewhere. dbSNP files can also be discarded or moved.

# Subset EUR and AFR samples to cut down on processing times downstream (move the output from each command to a separate folder). The files with all samples can be deleted or moved afterwards.
for f in *.vcf.gz;do bcftools view -S EUR_503_Both_Sexes -Oz -o EUR_503_Both.${f} $f --threads 10;done
for f in *.vcf.gz;do bcftools view -S AFR_660_Both_Sexes -Oz -o AFR_660_Both.${f} $f --threads 10;done

# Make a list of the vcf files for each superpopulation
ls *vcf.gz > vcf_list

# Manually reorganize each list so that it is ordered from chr1 to chr22 (autosomes)

# Concatenate SNVs on each chromosome into one file
bcftools concat --file-list vcf_list --output chr1-22.EUR_503_Both_rsID.phase3.genotypes.vcf.gz --threads 10 --naive
bcftools concat --file-list vcf_list --output chr1-22.AFR_660_Both_rsID.phase3.genotypes.vcf.gz --threads 10 --naive

# Move these concatenated files to the main 'SNV_Genotypes' directory. Everything else can be moved or deleted.


################### Filter analysis-grade SNVs with vcftools

# Go up to the main genotype directory
cd SNV_Genotypes

# Update vcf headers for compatibility with vcftools
sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' chr1-22.EUR_503_Both_rsID.phase3.genotypes.vcf > chr1-22.EUR_503_Both_rsID.4.2.phase3.genotypes.vcf
sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' chr1-22.AFR_660_Both_rsID.phase3.genotypes.vcf > chr1-22.AFR_660_Both_rsID.4.2.phase3.genotypes.vcf

# Filter EUR SNVs
vcftools --vcf chr1-22.EUR_503_Both_rsID.4.2.phase3.genotypes.vcf --out Filtered_EUR_BOTH_SEXES_chr1-22 --min-alleles 2 --max-alleles 2 --remove-indels --maf 0.01 --hwe .000001 --max-missing 1 --keep EUR_358_GEUV_and_Phase3_hg38 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --recode --recode-INFO-all

# Filter AFR SNVs
vcftools --vcf chr1-22.AFR_660_Both_rsID.4.2.phase3.genotypes.vcf --out Filtered_YRI_BOTH_SEXES_chr1-22 --min-alleles 2 --max-alleles 2 --remove-indels --maf 0.01 --hwe .000001 --max-missing 1 --keep YRI_86_GEUV_and_Phase3_hg38 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --recode --recode-INFO-all

# compress filtered SNVs
bgzip -@ 10 Filtered_EUR_BOTH_SEXES_chr1-22.recode.vcf
bgzip -@ 10 Filtered_YRI_BOTH_SEXES_chr1-22.recode.vcf


################### Genotype PCA Analysis with Plink

# Prune the filtered, EUR genotype data
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf Filtered_EUR_BOTH_SEXES_chr1-22.recode.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out Plink_EUR_PCA/Filtered_EUR_BOTH_SEXES_pruned

# Prune the filtered, YRI genotype data
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf Filtered_YRI_BOTH_SEXES_chr1-22.recode.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out Plink_YRI_PCA/Filtered_YRI_BOTH_SEXES_pruned

# Run PCA analysis with Plink, on pruned EUR data
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf Filtered_EUR_BOTH_SEXES_chr1-22.recode.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract Plink_EUR_PCA/Filtered_EUR_BOTH_SEXES_pruned.prune.in --pca var-wts --out Plink_EUR_PCA/Filtered_EUR_BOTH_SEXES_pca

# Run PCA analysis with Plink, on pruned YRI data
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf Filtered_YRI_BOTH_SEXES_chr1-22.recode.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract Plink_YRI_PCA/Filtered_YRI_BOTH_SEXES_pruned.prune.in --pca var-wts --out Plink_YRI_PCA/Filtered_YRI_BOTH_SEXES_pca


################### Convert VCFs to Plink format and generate 0/1/2 matrices for MatrixeQTL

# Make plink bed files using the filtered genotype data
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf Filtered_EUR_BOTH_SEXES_chr1-22.recode.vcf.gz --set-missing-var-ids @:# --keep-allele-order --make-bed --out Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf Filtered_YRI_BOTH_SEXES_chr1-22.recode.vcf.gz --set-missing-var-ids @:# --keep-allele-order --make-bed --out Plink_YRI_Genotypes/Filtered_YRI_BOTH_SEXES_chr1-22.recode.plink

# Make a list of duplicate RsIDs from the plink.bim files (will be used downstream for LD clumping top SNVs from eQTL analyses)
cut -f 2 Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink.bim | sort | uniq -d > Plink_EUR_Genotypes/duplicate_rsids.txt
cut -f 2 Plink_YRI_Genotypes/Filtered_YRI_BOTH_SEXES_chr1-22.recode.plink.bim | sort | uniq -d > Plink_YRI_Genotypes/duplicate_rsids.txt

# Make plink.raw files, which will serve as the base for the 0/1/2 matrices
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink --keep-allele-order --recodeA --out Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink.012
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile Plink_YRI_Genotypes/Filtered_YRI_BOTH_SEXES_chr1-22.recode.plink --keep-allele-order --recodeA --out Plink_YRI_Genotypes/Filtered_YRI_BOTH_SEXES_chr1-22.recode.plink.012

# From the newly generated .raw files, remove columns 2-6
gcut --complement -d ' ' -f 2-6 Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink.012.raw > Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink.012
gcut --complement -d ' ' -f 2-6 Plink_YRI_Genotypes/Filtered_YRI_BOTH_SEXES_chr1-22.recode.plink.012.raw > Plink_YRI_Genotypes/Filtered_YRI_BOTH_SEXES_chr1-22.recode.plink.012

# On the newly generated .012 files, remove the first row (RsIDs). ONLY RUN ONCE SINCE THESE COMMANDS ACT DIRECTLY ON THE FILE.
sed -i '' '1d' Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink.012
sed -i '' '1d' Plink_YRI_Genotypes/Filtered_YRI_BOTH_SEXES_chr1-22.recode.plink.012

# Run the 'SNV_Genotype_Preparation_Supplement.R' script, which will provide new SNP IDs to use (since RsIDs may be redundant)

# Append the new SNP IDs to the 0/1/2 matrices, using temporary files
cat Plink_EUR_Genotypes/INPUT_Dependency_SNP_IDs_row_EUR.txt Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink.012 > Plink_EUR_Genotypes/temp1
cat Plink_YRI_Genotypes/INPUT_Dependency_SNP_IDs_row_YRI.txt Plink_YRI_Genotypes/Filtered_YRI_BOTH_SEXES_chr1-22.recode.plink.012 > Plink_YRI_Genotypes/temp2

# Transpose the temporary files to generate the final 0/1/2 matrices in MatrixeQTL format (rows=SNP_IDs, columns=Samples)
datamash -W transpose < Plink_EUR_Genotypes/temp1 > Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink.012
datamash -W transpose < Plink_YRI_Genotypes/temp2 > Plink_YRI_Genotypes/Filtered_YRI_BOTH_SEXES_chr1-22.recode.plink.012
