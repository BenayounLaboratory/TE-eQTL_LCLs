README -- RNASeq Analysis of Candidate Regulator OE
############################

  1. Trim reads with fastp, check quality with FastQC, and map to hg38 with STAR (Trim_and_Map_GM12878_OE.sh)

  2. Count reads mapping to genes or repeats using TETranscripts (TETranscripts_GM12878_OE.sh)

  3. Filter lowly expressed genes, run DESeq2 for DE genes/TEs comparing empty vector vs overexpression vectors, and run MDS/PCA analyses (Process_RNASeq_GM12878_OE.R, Process_RNASeq_functions.R)

  4. Run GSEA enrichment analysis:
      - Run GSEA and plot the top 5 gene sets in each direction (Run_GSEA_GM12878_OE.R)
          - Note: The gene sets prepared for the GSEA analysis in '2_RNASeq_Analysis_Top_eQTLs' are used here
          - Note: The functions ('Run_GSEA_Against_Top_eQTLs_Functions.R') prepared for the GSEA analysis in '2_RNASeq_Analysis_Top_eQTLs' are used here
      - Find overlapping significant gene sets between OE conditions (Run_Comparison_of_GSEA_Results.R)
          - Note: The functions ('Run_GSEA_Against_Top_eQTLs_Functions.R') prepared for the GSEA analysis in '2_RNASeq_Analysis_Top_eQTLs' are used here
