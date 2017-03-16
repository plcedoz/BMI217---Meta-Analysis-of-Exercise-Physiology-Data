#This is a pilot study based on GSE28422 to explore the effect of training on gene expression

library(gplots)
library(dplyr)

rm(list = ls())

studies = c("GSE28422_1", "GSE28422_2", "GSE28422_3")

for (study in studies){
  assign(study, read.table(sprintf("gene_lists/%s", study), header = TRUE))
}

datasets = list()
for (study in studies){
  datasets[[study]] = read.table(sprintf("gene_lists/%s", study), header = TRUE)
  datasets[[study]] = filter(datasets[[study]], logFC>0.5)
}

split_multiple_gene_names <- function(geneL){
  for(i in 1: length(geneL)){
    if((regexpr('\\///',geneL[i])) != -1){
      cur_gene <- geneL[i]
      gene1 <- substr(cur_gene, 1, regexpr('\\///',cur_gene)-1)
      geneL = c(geneL, gene1)
      gene2 <- substr(cur_gene, regexpr('\\///',cur_gene) + 3, nchar(cur_gene))
      if ((regexpr('\\///',gene2)) != -1){
        cur_gene2 = gene2
        gene2 <- substr(cur_gene2, 1, regexpr('\\///',cur_gene2)-1)
        gene3 <- substr(cur_gene2, regexpr('\\///',cur_gene2)+3, nchar(cur_gene2))
        geneL = c(geneL, gene2)
        geneL = c(geneL, gene3)
      }
      else{
        geneL = c(geneL, gene2)
      }
      geneL = geneL[geneL != cur_gene]
    }
  }
  return(geneL)
}

gene_lists = list()
for (study in studies){
  gene_lists[[study]] = unique(as.character(datasets[[study]]$Gene.symbol))
  gene_lists[[study]] = split_multiple_gene_names(gene_lists[[study]])
}

all_genes = c()
for (study in studies){
  all_genes = union(all_genes, gene_lists[[study]])
}

heatmap_matrix = matrix(0, nrow = length(all_genes), ncol = length(studies))
rownames(heatmap_matrix) = all_genes
colnames(heatmap_matrix) = studies
for (study in studies){
  heatmap_matrix[,study] = all_genes %in% gene_lists[[study]]
}
colnames(heatmap_matrix) = c("acute", "chronic1", "chronic2")
heatmap.2(heatmap_matrix, trace = "none", col = (c("white", "black")), margins = c(4,4), cexCol = 1)
most_common_genes = cbind(heatmap_matrix,rowSums(heatmap_matrix))
most_common_genes = most_common_genes[order(most_common_genes[,ncol(most_common_genes)], decreasing = TRUE),]


#Create a data matrix containing the qvalue of each genes for each studies for a different heatmap
heatmap_matrix_2 = matrix(0, nrow = length(all_genes), ncol = length(studies))
rownames(heatmap_matrix_2) = all_genes
colnames(heatmap_matrix_2) = studies
for (study in studies){
  heatmap_matrix_2[,study] = all_genes %in% gene_lists[[study]]
  qvalue = datasets[[study]]$adj.P.Val
  index = 1
  for (i in 1:nrow(heatmap_matrix_2)){
    if (heatmap_matrix_2[i,study] == 1){
      heatmap_matrix_2[i,study] = qvalue[index]
      index = index + 1
    }
  } 
}
heatmap_matrix_2[heatmap_matrix_2 == 0] = 1
heatmap_matrix_2 = log10(heatmap_matrix_2)

#Plot the heatmaps
colnames(heatmap_matrix_2) = c("acute", "chronic1", "chronic2")
heatmap.2(heatmap_matrix_2, trace = "none", col = colorRampPalette(c("black","white"))(100), symbreaks = FALSE, margins = c(6,4), cexCol = 1.5)



