# dependencies:
# install.packages("ggplot2")
# install.packages("gridExtra")
# install.packages("plotly")
suppressPackageStartupMessages(library("plotly"))

# path to the DiffBind table with genes

# read the file into a dataframe
diff_df <- read.table('GSE58559.txt', 
                      fill=T, sep="\t", stringsAsFactors=F, header=T, row.names=1, 
                      colClasses=c('character'))   

# check some attributes of the data
colnames(diff_df)
dim(diff_df)
# keep only the fields needed for the plot
# FDR = false discovery rate = adjusted p value = significance 
# We can use P-Value rather than FDR
diff_df <- diff_df[c("Gene.symbol", "logFC", "P.Value")]

# preview the dataset; data required for the plot
head(diff_df)
# In order to color the points on the plot easier, we will add an extra column to 
# the data to mark whether each entry passes our cutoffs for Fold change and Significance. 
# Our criteria for this experiment were:
# Fold change > 1.5x change
# FDR < 0.5
# add a grouping column; default value is "not significant"
diff_df["group"] <- "NotSignificant"

# for our plot, we want to highlight 
# FDR < 0.05 (significance level)
# Fold Change > 1.5

diff_df$P.Value <- sapply(diff_df$P.Value, as.numeric)
diff_df$logFC <- sapply(diff_df$logFC, as.numeric)

# change the grouping for the entries with significance but not a large enough Fold change
diff_df[which(diff_df['P.Value'] < 0.05 & abs(diff_df['logFC']) < 1.5 ),"group"] <- "Significant"

# change the grouping for the entries a large enough Fold change but not a low enough p value
diff_df[which(diff_df['P.Value'] > 0.05 & abs(diff_df['logFC']) > 1.5 ),"group"] <- "FoldChange"

# change the grouping for the entries with both significance and large enough fold change
diff_df[which(diff_df['P.Value'] < 0.05 & abs(diff_df['logFC']) > 1.5 ),"group"] <- "Significant&FoldChange"
# Find and label the top peaks..
top_peaks <- diff_df[with(diff_df, order(logFC, P.Value)),][1:5,]
top_peaks <- rbind(top_peaks, diff_df[with(diff_df, order(-logFC, P.Value)),][1:5,])


# Add gene labels for all of the top genes we found
# here we are creating an empty list, and filling it with entries for each row in the dataframe
# each list entry is another list with named items that will be used by Plot.ly
a <- list()
for (i in seq_len(nrow(top_peaks))) {
  m <- top_peaks[i, ]
  a[[i]] <- list(
    x = m[["logFC"]],
    y = -log10(m[["P.Value"]]),
    text = m[["Gene.symbol"]],
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 0.5,
    ax = 20,
    ay = -40
  )
}

# make the Plot.ly plot
p <- plot_ly(data = diff_df, x = diff_df$logFC, y = -log10(diff_df$P.Value), 
             text = diff_df$Gene.symbol, mode = "markers", color = diff_df$group) %>% 
  layout(title ="Volcano Plot") %>%
  layout(annotations = a)
p
