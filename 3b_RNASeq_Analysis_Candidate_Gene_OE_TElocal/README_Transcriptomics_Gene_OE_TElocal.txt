README -- RNASeq Analysis of Candidate Regulator OE (TElocal)
############################

  1. Count reads mapping to genes or repeats using TElocal (TElocal_GM12878_OE.sh)

  2. Filter lowly expressed genes, run DESeq2 for DE genes/TEs comparing empty vector vs overexpression vectors, and run MDS/PCA analyses (Process_RNASeq_GM12878_OE.R)

  3. Run GSEA enrichment analysis:
      - Run GSEA using the stratified repeat gene sets (Run_GSEA_GM12878_OE.R)
      - Generate consolidated bubble plots for the stratified repeat GSEA results (Compare_repeats_for_each_treatment.R)

