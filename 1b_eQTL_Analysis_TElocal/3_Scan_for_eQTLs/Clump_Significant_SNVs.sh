# Go to directory with the significant SNVs
cd Significant_SNVs

# Clump distal intergenic L1 trans-eQTL SNVs using Plink
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile /Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink --exclude /Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/duplicate_rsids.txt --clump L1_intergenic_distal_trans_eQTLs_SIGNIFICANT_EUR_2023-11-15.txt --clump-p1 0.05 --clump-p2 1 --clump-r2 0.10 --clump-kb 500 --clump-snp-field RsID --clump-field pvalue --set-missing-var-ids @:# --out L1_intergenic_distal_trans_eQTLs_SIGNIFICANT_EUR

# Clump nearby intergenic L1 trans-eQTL SNVs using Plink
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile /Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink --exclude /Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/duplicate_rsids.txt --clump L1_intergenic_nearby_trans_eQTLs_SIGNIFICANT_EUR_2023-11-15.txt --clump-p1 0.05 --clump-p2 1 --clump-r2 0.10 --clump-kb 500 --clump-snp-field RsID --clump-field pvalue --set-missing-var-ids @:# --out L1_intergenic_nearby_trans_eQTLs_SIGNIFICANT_EUR

# Clump intronic L1 trans-eQTL SNVs using Plink
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile /Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink --exclude /Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/duplicate_rsids.txt --clump L1_intronic_trans_eQTLs_SIGNIFICANT_EUR_2023-11-15.txt --clump-p1 0.05 --clump-p2 1 --clump-r2 0.10 --clump-kb 500 --clump-snp-field RsID --clump-field pvalue --set-missing-var-ids @:# --out L1_intronic_trans_eQTLs_SIGNIFICANT_EUR

# Clump exonic L1 trans-eQTL SNVs using Plink
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile /Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/Filtered_EUR_BOTH_SEXES_chr1-22.recode.plink --exclude /Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_EUR_Genotypes/duplicate_rsids.txt --clump L1_exonic_trans_eQTLs_SIGNIFICANT_EUR_2023-11-15.txt --clump-p1 0.05 --clump-p2 1 --clump-r2 0.10 --clump-kb 500 --clump-snp-field RsID --clump-field pvalue --set-missing-var-ids @:# --out L1_exonic_trans_eQTLs_SIGNIFICANT_EUR
