# 02_Interpret_heatmap.R
# Auth: james.opzoomer@kcl.ac.uk
# Desc: Take clusters from heatmap and manually metacluster
# Date: 12/02/2018

# Need to trace modify the write.FCS fucntion to alter the $PnR - done
# line 35 -> change 1025
## Helper functions
read_dfToFcs <- function(dat) {
  
  dta = dat
  
  head(dta)
  
  # you need to prepare some metadata
  meta = data.frame(name=dimnames(dta)[[2]],
                     desc=paste('this is column',dimnames(dta)[[2]],'from your CSV')
  )
  meta$range = apply(apply(dta,2,range),2,diff)
  meta$minRange = apply(dta,2,min)
  meta$maxRange = apply(dta,2,max)
  
  head(meta)
  
  mat.dta = as.matrix(dta)
  
  # all these are required for the following steps to work
  # a flowFrame is the internal representation of a FCS file
  ff = new("flowFrame",
            exprs=mat.dta,
            parameters=AnnotatedDataFrame(meta))
  
  out <- paste0(gswd, 'FCS_clusterFile', '_', pid)
  
  # now you can save it back to the filesystem
  # Need to find a permanent way to edit the filesystem
  write.FCS(ff,paste(out,'FCS',sep='.'))
}


## End of helper functions

library(flowCore)
library(ggplot2)

# Load config file
source('~/Desktop/Arnold/Data/2018_01_Louvain_clustering/FCS_clustering/00_panel_config.R')

setwd(gswd)

# LOOK AT THE HEATMAP AND FIGURE OUT WHAT YOU WHAT LEVEL YOU WANT TO CUT THE HEATMAP AT 
# variable to cut the heatmap
hCut = 4

# Load the expression matricx with cluster assignment
n = paste0(gswd, 'scaled_avg_expression')
load(n)

# New heatmap object with metacluster grouping
map = pheatmap(scaled_avg_expression, cutree_rows = hCut, cellheight = 8, cellwidth = 8, fontsize_col = 8, fontsize_row = 8, filename = "heatmap.png")

# Save the hclust object
n = paste0(gswd, 'hClust_object')
save(map, file = n)
# load(n)

# separate the hclust objects
# create refrence list of the metaclusters and their children
map.clust <- cbind(scaled_avg_expression, cluster = cutree(map$tree_row, k = hCut))
map.clust = as.data.frame(map.clust)

# Generate the metacluster column
total_clust_MSC$metaCluster =  as.factor(total_clust_MSC$phenograph_cluster)

# Map values
total_clust_MSC$metaCluster = mapvalues(total_clust_MSC$metaCluster, from = rownames(map.clust), to = map.clust$cluster)

# Save this matrix
n = 'cell.matrix.with.cell.cluster.assignment'
save(total_clust_MSC, file = n)

# Vizualise, note that meta_pheno_cluster is switched
ggplot(total_clust_MSC, aes(x=tSNE_x, y=tSNE_y, col=metaCluster)) + geom_point(size = 0.6) + theme_bw()

fp = paste0(gswd, "MSC_metacluster.png")
ggsave(fp)

# Dope Viz
ggplot(total_clust_MSC, aes(x=tSNE_x, y=tSNE_y, col=CD206)) + geom_point(size = 0.6) + theme_bw() + scale_color_gradientn("CD206", colours = colorRampPalette((rev(brewer.pal(n = 11, name = "Spectral"))))(50))
ggplot(total_clust_MSC, aes(x=tSNE_x, y=tSNE_y, col=MHCII)) + geom_point(size = 0.6) + theme_bw() + scale_color_gradientn("MHCII", colours = colorRampPalette((rev(brewer.pal(n = 11, name = "Spectral"))))(50))

# Produce the output report of cluster data
dfSummaryData = as.data.frame(table(total_clust_MSC$metaCluster))
colnames(dfSummaryData) = c('meta_cluster', 'count')
dfSummaryData$percent_total = (dfSummaryData$count / sum(dfSummaryData$count))*100
n = paste0(gswd, pid, '_summaryStats')
save(dfSummaryData, file = n)
write.csv(dfSummaryData, paste0(n, '.csv'))

# Convert this file back into an FCS file
# Need to convert the original non-transformed data
# Load original data
dfFCS_data = read_csv(fcs_file)
colnames(dfFCS_data) = dfPanel$markers
dfFCS_data$metaCluster = as.integer(total_clust_MSC$metaCluster)

dfInt_FCS = sapply(dfFCS_data, function(x){as.integer(x)})
read_dfToFcs(dfInt_FCS)



