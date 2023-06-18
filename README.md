# Double_Mutant_Gene_Expression_Analysis
Scripts used in Napier-Jameson et al., 2023 to identify transcriptomic changes that are unique to double mutants (more extreme than expected based on single mutant transcriptomes)

First, outputs of DESeq (gene expression analysis) or JUM (alternative splicing analysis) are processed through either the gene_expression_script.R (for gene expression) or splicing_script.R (for splicing analysis)

Second, these outputs are run through unique_dysregulation_identification.R to find genes or splicing events that are "uniquely dysregulated" (i.e. more extreme than predicted by an additive model based on the two single mutants)
