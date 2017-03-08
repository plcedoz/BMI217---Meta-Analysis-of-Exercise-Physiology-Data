library(gplots)
library(dplyr)

studies = list.files(path = "gene_lists")

datasets = list()
for (study in studies){
  assign(study, read.table(sprintf("gene_lists/%s", study), header = TRUE))
}

for (study in studies){
  datasets[[study]] = read.table(sprintf("gene_lists/%s", study), header = TRUE)
  datasets[[study]] = filter(datasets[[study]], logFC > 1.3)
}

gene_lists = list()
for (study in studies){
  gene_lists[[study]] = unique(as.character(datasets[[study]]$Gene.symbol))
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



