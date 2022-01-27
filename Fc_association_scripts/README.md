Script Guide:

- Imputation and intial QC of data saved separately for now

- QC_PCA.R script contains: variable recoding for PLINK, QC thresholds, PCA analysis
- duplicate sorting.R contains code for recoding sample IDs, VAST and PATCH had the same IDs (was quite the rigamarol)
- association_testing.R shows some intial regression models run in PLINK/R as well as my attempts to input the results into FUMA for LD pruning etc
- pos2rsID.R is a script for converting the vcf coordinates back into rsIDs of GChr38
- snp_annotating.R contains a script for downloading uscb genome track browser data, and subsequent LD pruning attempts based on it
- Fc_associations.R contains assoc testing of Fcy region and attempts to plot the initial most significant snps by case control and genotype

