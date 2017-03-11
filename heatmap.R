library(gplots)
library(dplyr)

rm(list = ls())

studies = list.files(path = "gene_lists")
studies = studies [-16]
studies = studies [studies != "GSE28422.txt"]
studies = studies [studies != "GSE34788_2.txt"]
studies = studies [studies != "GSE58559.txt"]
studies = studies [studies != "GSE18583.txt"]
studies = studies [studies != "GSE60655.txt"]
studies = studies [studies != "GSE9405_1.txt"]
studies = studies [studies != "GSE9405_3.txt"]


for (study in studies){
  assign(study, read.table(sprintf("gene_lists/%s", study), header = TRUE))
}

datasets = list()
for (study in studies){
  datasets[[study]] = read.table(sprintf("gene_lists/%s", study), header = TRUE)
  datasets[[study]] = filter(datasets[[study]], logFC > 0.5)
}

split_multiple_gene_names <- function(geneL){
  for(i in 1: length(geneL)){
    print(geneL[i])
  }
}

gene_lists = list()
for (study in studies){
  gene_lists[[study]] = unique(as.character(datasets[[study]]$Gene.symbol))
  split_multiple_gene_names(gene_lists[[study]])
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

heatmap.2(heatmap_matrix, trace = "none", col = (c("white", "black")))
most_common_genes = cbind(heatmap_matrix,rowSums(heatmap_matrix))
most_common_genes = most_common_genes[order(most_common_genes[,11], decreasing = TRUE),]


#Clustering: is it possible to find acute vs chronic?
heatmap_matrix = t(heatmap_matrix)
dataframes_clusters = kmeans(heatmap_matrix, centers = 2, nstart = 10)
dataframes_clusters


