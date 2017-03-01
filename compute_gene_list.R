library(GEOquery)
library(limma)
library(dplyr)

rm(list =ls())

source("get_data_gse.R")

gse_ref = "GSE8668"
data_gse = get_data_gse(gse_ref)
data_matrix = data_gse$data_matrix
gsm_description = data_gse$gsm_description

#We build the design matrix (case/control) from the gsm_description file
after_exercise = gsm_description$V3 == "after-exercise"
design = model.matrix(~ after_exercise)
rownames(design) = gsm_description$V1
fit <- lmFit (data_matrix, design)
fit <- eBayes (fit)
#Compute genes are differentially expressed at an adjusted p-value of ???? If we choose 0.05: 0 genes
toptable <- toptable(fit, coef = 2,  p.value = 1, number = 50000)
list_of_genes = select(toptable, ID, P.Value, adj.P.Val)

write.table(list_of_genes[1:1000,], sep = "\t", file = "gene_list_GSE8668", quote = FALSE, row.names = FALSE)
#If we use Fold Change we don't have any genes!
#toptable_lfc <- toptable(fit, coef = 2,  p.value = 0.5, lfc = 2^1.3, number = 50000)
#number_lfc = dim(toptable_lfc)[1]



