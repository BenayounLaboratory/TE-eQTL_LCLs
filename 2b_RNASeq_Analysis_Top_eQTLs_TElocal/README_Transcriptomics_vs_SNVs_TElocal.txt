README -- RNASeq Analysis of Top eQTLs (TElocal)
############################

  1. Run DESeq2 analysis for DE genes/TEs as a function of the most significant L1 trans-eQTLs in each clump identified in '1_eQTL_Analysis_TElocal' (Run_Genotype_DESeq.R)

  2. Run GSEA enrichment analysis:
      - Prepare stratified repeat GMT gene sets for GSEA (Prepare_Gene_Sets.R)
      - Run GSEA using the stratified repeat gene sets (Run_GSEA_Against_Top_eQTLs.R)
      - Run pathway GSEA for the ZSCAN26 SNV and plot the top 5 gene sets in each direction (Run_GSEA_TElocal_New_eQTL_Pathways.R)
      - Generate consolidated bubble plots for the stratified repeat GSEA results (Compare_repeats_for_each_SNV.R)
     