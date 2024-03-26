README -- eQTL analysis
############################

  1. Prepare SNV (single nucleotide variant) and SV (structural variant) genotype data for eQTL and transcriptomic analyses
      - Filter SNV genotypes, convert VCF files to Plink bed and 0/1/2 formats, and prune SNV genotypes for PCA analysis (SNV_Genotype_Preparation_Instructions.sh, SNV_Genotype_Preparation_Supplement.R)
      - Filter SV genotypes and prune SV genotypes for PCA analysis (SV_Genotype_Preparation_Instructions.sh)
      - Quantify L1 and Alu insertions and deletions (TE_Copy_Number_Quantification_Instructions.sh, TE_Copy_Number_Quantification_Supplement.R)
      - Plot SNV genotype PCA results and export PCs as a text file (Run_Population_Structure_Analysis_SNVs.R)
      - Plot SV genotype PCA results and export PCs as a text file (Run_Population_Structure_Analysis_SVs.R)

  2. Trim reads with fastp, check quality with FastQC, and map to hg38 with STAR (Trim_and_Map_Geuvadis.sh)

  3. Count reads mapping to genes or repeats using TETranscripts (TETranscripts_Geuvadis.sh)

  4. Prepare RNAseq data for eQTL analysis
      - Filter lowly expressed genes, VST transform counts, remove batch effects, and inverse normal transform the expression values (Process_RNASeq_for_eQTL.R, Process_RNASeq_for_eQTL_functions.R)
      - As a parallel/supplemental approach, take VST expression data, identify PEER factors, remove known batch effects and PEER factors, and inverse normal transform the expression values (Process_RNASeq_for_eQTL_PEER.R, Process_RNASeq_for_eQTL_functions.R)
      - Plot known batches (Plot_Batch_Effects.R, Plot_Batch_Effects_Functions.R)

  5. Define significant eQTLs
      - Prepare additional files needed for eQTL analysis, including a gene position file and a sample permutations file (Prepare_Remaining_Input_Files.R, Prepare_Remaining_Input_Files_Functions.R)
      - Run gene cis-eQTL and L1 trans-eQTL analyses, including analyses with sample permutations (Run_MatrixeQTL_v2.R, Run_MatrixeQTL_Functions.R)
      - Calculate empirical FDR thresholds using permutation results, and extract significant eQTLs (Define_thresholds_and_extract_significant_SNVs_v2.R, Define_thresholds_and_extract_significant_SNVs_functions.R)
      - Clump eQTLs by p-value (Clump_Significant_SNVs_v2.sh)

  6. Combine cis- and trans- eQTL results
      - Combine cis/trans eQTLs sharing an SNV, run linear regressions on linked genes/TEs, and generate box whisker plots for unique SNV-gene-TE trios (Integrate_eQTLs_and_Generate_Plots_v2.R, Integrate_eQTLs_and_Generate_Plots_Functions.R)
      - Generate annotated Manhattan plots (Generate_Manhattan_Plots_v2.R)
      - Plot EBV vs IL16 correlations pre- and post-correction of eQTL data (Check_EBV_Correlations.R)

  7. Evaluate top EUR cis/trans eQTLs in an independent YRI population, extract significant eQTLs, combine cis/trans eQTLs, and generate box-whisker plots (Check_Replication_of_eQTLs.R, Check_Replication_of_eQTLs_Functions.R, Clump_Significant_SNVs.sh)

  8. Evaluate whether cis-genes mediate SNV-TE relationships (Run_Mediation_Analysis.R)

  9. Evaluate whether there is an enrichment of TE copies near significant, index eQTL SNVs compared to random (Run_check_TE_Enrichments.R)
