library(gplots)
library(dplyr)

rm(list = ls())

studies = list.files(path = "gene_lists")
studies = studies [studies != "gene_list_description"]

studies = studies [studies != "GSE11761.txt"]
studies = studies [studies != "GSE44051.txt"]
studies = studies [studies != "GSE8668.txt"]
studies = studies [studies != "GSE28422_1"]
studies = studies [studies != "GSE58559_1.txt"]
studies = studies [studies != "GSE58559_2.txt"]
studies = studies [studies != "GSE58559_3.txt"]
studies = studies [studies != "GSE68072.txt"]
studies = studies [studies != "GSE59088_1.txt"]
studies = studies [studies != "GSE59088_2.txt"]


#for (study in studies){
#  assign(study, read.table(sprintf("gene_lists/%s", study), header = TRUE))
#}

datasets = list()
for (study in studies){
  datasets[[study]] = read.table(sprintf("gene_lists/%s", study), header = TRUE)
  datasets[[study]] = filter(datasets[[study]], adj.P.Val < 0.05)
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
  #gene_lists[[study]] = split_multiple_gene_names(gene_lists[[study]])
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
heatmap.2(heatmap_matrix_2, trace = "none", col = colorRampPalette(c("black","white"))(100), symbreaks = FALSE)

write.csv(heatmap_matrix_2, "heatmap_matrix")

heatmap.2(heatmap_matrix, trace = "none", col = (c("white", "black")))
most_common_genes = cbind(heatmap_matrix,rowSums(heatmap_matrix))
most_common_genes = most_common_genes[order(most_common_genes[,ncol(most_common_genes)], decreasing = TRUE),]

final_gene_list = rownames(most_common_genes[most_common_genes[,ncol(most_common_genes)]>0,])

#Clustering: is it possible to find acute vs chronic?
heatmap_matrix_2 = t(heatmap_matrix_2)
dataframes_clusters = kmeans(heatmap_matrix_2, centers = 3, nstart = 10)
dataframes_clusters

