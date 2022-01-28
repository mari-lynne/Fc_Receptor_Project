#Thoughts ####

#Route 1
#Prune for LD
#Take our most significant snps
#Then look to see what they are in LD with in terms of SNPs with functional consequences
#Include these in the analysis


#Route2
#Take sig SNPs and include functional snps that are in LD
#Prune but prioritise, sig and functional snps

#Route 3
#Group SNPs that are in LD together and then do haplotype analysis 
#I think do this anyways

#check merging of strands/alleles?
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6099125/
In Illumina annotation each SNP is defined with design allele nucleotides, and these occur on the same strand as the probe sequence; the order in which the alternative alleles are given specifies the generic A and B allele designations [1]. To illustrate, for a SNP defined as [T/G], the A allele is T and the B allele is G. In Affymetrix allele-specific hybridization technology, the letter codes A and B are assigned differently and could therefore occur on either the probe or target strand

Im pretty sure I aligned the different datasets to the refence genome so strand using will rayner topmed perl script, so it should not be an issue, and it was also checked following the output of imputation against the fa of hg38 (takes hours) but I want to double check this especially considering the different platforms

#Do strand checker for each study seperately pre-imputation (v sure I did this but will double check) - used will rayners tool
#After strand checking see how things merge pre-imputation
