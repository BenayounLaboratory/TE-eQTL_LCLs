README -- eQTL analysis (TElocal)
############################

  0. Prepare GTF files for repeats stratified by genomic context
      - Convert Hammell lab repeat location files to GTF format (TElocal_locations_to_GTF_v1_Nov_21_2023.R)
      - Split the Hammell lab TE locus GTF file by genomic context

  1. Count reads mapping to genes or repeat loci using TElocal (TElocal_Geuvadis.sh)

  2. Prepare RNAseq data for eQTL analysis
      - Filter lowly expressed genes, VST transform counts, remove batch effects, and inverse normal transform the expression values (Process_RNASeq_for_eQTL.R)

  3. Define significant eQTLs
      - Run gene cis-eQTL and stratified L1 trans-eQTL analyses (Run_MatrixeQTL.R)
      - Extract significant eQTLs using BH FDR (Define_thresholds_and_extract_significant_SNVs.R)
      - Clump eQTLs by p-value (Clump_Significant_SNVs.sh)

  4. Combine cis- and trans- eQTL results
      - Combine cis/trans eQTLs sharing an SNV, run linear regressions on linked genes/TEs, and generate box whisker plots for unique SNV-gene-TE trios (Integrate_eQTLs_and_Generate_Plots.R)
      - Generate annotated Manhattan plots (Generate_Manhattan_Plots.R)

  5. Evaluate top EUR cis/trans eQTLs in an independent YRI population, extract significant eQTLs, combine cis/trans eQTLs, and generate box-whisker plots (Check_Replication_of_eQTLs.R, Clump_Significant_SNVs.sh)

  6. Evaluate whether cis-genes mediate SNV-TE relationships (Run_Mediation_Analysis.R)
