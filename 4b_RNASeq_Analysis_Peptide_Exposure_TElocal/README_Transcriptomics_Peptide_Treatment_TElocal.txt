README -- RNASeq Analysis of rhIL16 Peptide Exposure (TElocal)
############################

  1. Count reads mapping to genes or repeats using TElocal (TElocal_rhIL16.sh)

  2. Filter lowly expressed genes, run DESeq2 for DE genes/TEs comparing 0 hr vs 24/48 hr exposures, and run MDS/PCA analyses (Process_RNASeq_GM12878_rIL16_TElocal.R)

  3. Run GSEA enrichment analysis:
      - Run GSEA using the stratified repeat gene sets (Run_GSEA_GM12878_rIL16.R)
      - Generate consolidated bubble plots for GSEA results with stratified repeat gene sets (Compare_repeats_for_each_treatment.R)
