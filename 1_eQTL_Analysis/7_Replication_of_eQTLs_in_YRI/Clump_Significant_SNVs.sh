# Go to directory with the significant SNVs
cd Scan_for_eQTLs

# Clump L1 trans-eQTL SNVs using Plink
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile /Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_YRI_Genotypes/Filtered_YRI_BOTH_SEXES_chr1-22.recode.plink --exclude /Users/juanb/Desktop/2023_TE_eQTL_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/SNV_Genotypes/Plink_YRI_Genotypes/duplicate_rsids.txt --clump L1_subfamily_trans_eQTLs_YRI_FDR5_2023-05-23.txt --clump-p1 0.05 --clump-p2 1 --clump-r2 0.10 --clump-kb 500 --clump-snp-field RsID --clump-field pvalue --set-missing-var-ids @:# --out L1_subfamily_trans_eQTLs_YRI_FDR5
