README -- RNASeq Analysis of Top eQTLs
############################

  1. Run DESeq2 analysis for DE genes/TEs as a function of the most significant L1 trans-eQTLs in each clump identified in '1_eQTL_Analysis' (Run_Genotype_DESeq.R, Run_Genotype_DESeq_functions.R)

  2. Run GSEA enrichment analysis:
      - Prepare GMT gene sets for GSEA (Prepare_Gene_Sets.R)
      - Run GSEA and plot the top 5 gene sets in each direction (Run_GSEA_Against_Top_eQTLs.R, Run_GSEA_Against_Top_eQTLs_Functions.R)
      - Find overlapping significant gene sets changing with the top 3 eQTLs with a cis-mediator (Compare_Top_3_Independent_SNVs.R)
      - Find overlapping significant gene sets changing with the top 2 eQTLs without a cis-mediator (Compare_Top_2_SNVs_without_mediator.R)
      - Find overlapping significant gene sets changing with all top 5 eQTLs (Compare_Top_5_Independent_SNVs.R)
