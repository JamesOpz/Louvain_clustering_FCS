# 01_FCS_clustering.R
# Auth: james.opzoomer@kcl.ac.uk
# Desc: tSNE dimensionality reduction followed by Louvain clustering
# Date: 29/01/2018

# Load config file
source('~/Desktop/Arnold/Data/2018_01_Louvain_clustering/FCS_clustering/00_panel_config.R')

# TO DO:
# I want to be able to do this with mutiple files and autogate quantify
# Neet to sort columns to assign markers
# Need to match markers (columns to subset data.frame)

list.of.packages <- c("ggplot2", "Rphenograph","flowCore","pheatmap","Rtsne","plyr","Biobase","readr","dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load dependencies
library(Rphenograph)
library(flowCore)
library(ggplot2)
library(pheatmap)
library(Rtsne)
library(plyr)
library(Biobase)
library(readr)
library(dplyr)

setwd(gswd)

# Load the FCS files
fcs_file
fp = paste0(gswd, fcs_file)
dfFCS_data = read_csv(fp)
colnames(dfFCS_data) = dfPanel$markers


# Scale the channel values by factor 150
dfTransform <- dfFCS_data[,3:n_markers]
dfTransform_asinh <- asinh(dfTransform*0.0066)


# re-alias
dfFCS_trans <- dfFCS_data
dfFCS_trans[,3:n_markers] <- dfTransform_asinh

# Convert to matrix
matMSC_trans <- as.matrix(dfFCS_trans)

# TD
# Subsample the FCS file
# Check how many FCS evensts there are and reduce if nessescary
set.seed(40)


# Need to work out the subsetting aspect of this.
# reset random seed for reporducibility
set.seed(420)
length(dfFCS_trans$FSC_A)
start.time <- Sys.time()
tsne_out <- Rtsne(dfFCS_trans[,3:n_markers], perpelexity = 65, theta = 0.3, pca = FALSE, max_iter = 1000)
end.time <- Sys.time()
time.taken <- end.time - start.time

time.taken

plots <- tsne_out$Y
colnames(plots) <- c("tSNE_x", "tSNE_y")
plots_df <- as.data.frame(plots)


dfFCS_trans$tSNE_x <- plots_df$tSNE_x
dfFCS_trans$tSNE_y <- plots_df$tSNE_y

# Viz the output of the tSNE
ggplot(dfFCS_trans, aes(x= tSNE_x, y= tSNE_y))+ 
  geom_point(shape = ".") + 
  theme_bw()

fp = paste0(gswd, "tSNE_raw.png")
ggsave(fp)

dfFCS_trans_mat = as.matrix(dfFCS_trans)
Rphenograph_out <- Rphenograph(dfFCS_trans_mat[,3:n_markers], k = 50)

dfFCS_trans$phenograph_cluster <- factor(membership(Rphenograph_out[[2]]))

# Viz the phenograph results
ggplot(dfFCS_trans, aes(x=tSNE_x, y=tSNE_y, col=phenograph_cluster)) + geom_point(size = 0.6) + theme_bw()

fp = paste0(gswd, "tSNE_pheno_clusters.png")
ggsave(fp)

total_clust_MSC = dfFCS_trans

n = paste0(gswd, 'clusterMatrix')
save(total_clust_MSC, file = n)
# Load statement for stopping point
# load(n)

# Make the phenograph average protein expression for each marker
# Need to make generic based off markers selected

scale_MSC = lapply(total_clust_MSC[,3:n_markers], function(x){ 
  
  normalized = (x-quantile(x, 0.01))/(quantile(x, 0.99)-quantile(x, 0.01))
})

scale_df <- as.data.frame(scale_MSC)

scale_df$phenograph_cluster = total_clust_MSC$phenograph_cluster  

## Without scaling
# Avg_expression = total_clust_MSC %>% group_by(phenograph_cluster) %>% summarise_all(funs(median))

## With scaling
Avg_expression = scale_df %>% group_by(phenograph_cluster) %>% summarise_all(funs(median))

scaled_avg_expression =  Avg_expression[,markers[3:n_markers]]

# Need to re-adapt the heatmap to scale the colours
# This should be a 0-1 scaling system

# scaled_avg_expression <- scale(Average_expression)

# Generate the heatmap with without trimming
# Inspect the heatmap and decide where to draw the line on your metaclusters
uncut_map = pheatmap(scaled_avg_expression, cellheight = 8, cellwidth = 8, fontsize_col = 8, fontsize_row = 8, filename = "noCutHeatmap.png")

n = paste0(gswd, 'scaled_avg_expression')
save(scaled_avg_expression, file = n)









