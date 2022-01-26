
# snps and gene expression T1T2

setwd("C:/Users/mjohnson/Documents/RNA")
Gene <- read.csv(file = "genotype_snps.csv")

sample_key <- read.csv(file = "sample_key_t1t2.csv")

Gene <- Gene %>% left_join(sample_key, by = "IID")

library(dplyr)

str(Gene)
Gene2 <- Gene %>% select(LabID, OUTCOME, GENOTYPE_462, GENOTYPE_11)

write.csv(Gene2, "Gene2.csv")

#counting genotype frequencies

GG <- Gene %>% filter(GENOTYPE_462 == "A G")

#recode genotypes

Gene2$GENOTYPE_462 <- ifelse(Gene2$GENOTYPE_462  == "A A", 0, ifelse(Gene2$GENOTYPE_462 == "A G", 1, ifelse(Gene2$GENOTYPE_462 == "G A", 1, ifelse(Gene2$GENOTYPE_462 == "G G", 2, NA))))

Gene2$GENOTYPE_11 <- ifelse(Gene2$GENOTYPE_11  == "G G", 0, ifelse(Gene2$GENOTYPE_11 == "A G", 1, ifelse(Gene2$GENOTYPE_11 == "G A", 1, ifelse(Gene2$GENOTYPE_11 == "A A", 2, NA))))

#restructure data frame
Gene3 <- t(Gene2)

#remove top header

write.csv(Gene3, "clean_table.csv")

#reload

clean_genes <- read.csv(file = "clean_table.csv")
#need to remove the xs at some point :( replace with P to match

library(dplyr)
library(stringr)

clean_genes <- clean_genes %>% rename_at(vars(starts_with("X")), 
          funs(str_replace(., "X", "P")))


### Expression Data ###

expression <- t(FCG_12FC)

write.csv(expression, "expression.csv")

full_list <- read.csv("Amber_expression.csv")

#search for genes and combine 

FCGR_expression <-filter(full_list, grepl('FCG', id))
FCRL_expression <-filter(full_list, grepl('FCRL', id))
FCAR_expression <-filter(full_list, grepl('FCAR', id))

FC_expression <- FCGR_expression %>% bind_rows(FCGR_expression, FCRL_expression,FCAR_expression)


##### Gene locations ###

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh = 37)
genes <- getBM(attributes=c('hgnc_symbol','chromosome_name','start_position','end_position'),mart = ensembl)

genes <- genes %>% rename(id = hgnc_symbol)

filtered_genes <- left_join(FC_expression, genes, by = "id")

filtered_genes <- filtered_genes %>% filter(chromosome_name == "1")

filtered_genes <- filtered_genes %>% dplyr::select(id, chromosome_name, start_position, end_position)

write.table(filtered_genes, "gene_loc.txt")
    

### snp locations ###

snp_loc <- read.csv(file = "snp_loc.csv")

write.table(snp_loc, "snp_loc.txt")

## match clean_genes to order of expression file 
test<- bind_rows(FC_expression, clean_genes) 

test<- FC_expression[names(clean_genes)] 

clean_genes %>% dplyr::all_of(names)

names <- colnames(FC_expression)
names2 <-colnames(clean_genes)

required_df <- clean_genes[clean_genes$names2 %in% FC_expression$names,]

required_df <-  FC_expression[FC_expression$names %in% clean_genes$names2,]

write.csv(FC_expression, file = "FC_expression.csv")

genes_matched <- read.csv("g.csv", fileEncoding="UTF-8-BOM")

str(genes_matched)

genes_matched <- genes_matched %>% dplyr::select(-LabID)

genes_matched <- genes_matched %>% rename(LabID = LabID.1)

### covariates ####

