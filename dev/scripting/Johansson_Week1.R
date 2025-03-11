

# Week 1 Homework - Introductory Analysis of the Data


# install.packages('Matrix')
library(Matrix)

# Loading helper functions
source('dev/utilities/dataprocessingHelperFunctions.R')

# if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

# Loading the flower data
a = load("data/flowerdata.RData", verbose=TRUE)

# Loading the original AraBOXcis network trained on bulk RNA-seq in seedlings
araboxcis <- read.csv(file = 'data/gboxNetwork22C.csv',header = TRUE)



# Contents of flower data
print(a)

# Dimensions of gbox
dim(gbox)

# Length of clust
length(clust)

# Names of the first 10 rows of gbox
rownames(gbox)[1:10]

# gbox column names
colnames(gbox)[1:10]

# First 50 rows and first 3 cols as matrix
as.matrix(gbox[1:50, 1:3])



# Cluster designations of the first 8 cells
clust[1:8]

# Cell counts in each cluster
table(clust)

# Plotting cell count against the clusters
plot(table(sort(clust)), 
     xlab = 'Cluster Name', 
     ylab = 'Number of Cells', 
     main = 'Histogram for the Flower Data')



# Dimensions of the Araboxcis matrix
dim(araboxcis)

# First 4 rows of the Araboxcis matrix
araboxcis[1:4,]

# Plotting a histogram of the edge scores
hist(araboxcis[,3])



# A list of the Transcription Factors 
tfs = unique(araboxcis[,1])

# Filtered list for only TFs also included in the single-cell RNA-seq dataset
tfSubs = tfs[which(tfs %in% rownames(gbox))]

# Number of TFs
length(tfSubs)



# Filter Genes and Cells With Very Low Values 

# Dimensions of the gbox
dim(gbox)

# Filtering-out cells with <1% of the genes expressed at all
thresh = 0.01
numberGenesPerCell = apply(gbox, 
                           2, 
                           function(i) {
                             length(which(i>0))
                             }
                           )
includeCells = which(numberGenesPerCell > (thresh * dim(gbox)[1])) 

gbox_filtered = gbox[,includeCells] 

dim(gbox_filtered)

# Removing genes that are expressed in <1% of the cells
numberCellPerGene = apply(gbox, 
                          1, 
                          function(i) {
                            length(which(i>0))
                            }
                          )
includeGenes = which(numberCellPerGene > (thresh * dim(gbox)[2]))

gbox_filtered = gbox_filtered[includeGenes,]

dim(gbox_filtered)



# First Visualisation With UMAP

install.packages('umap')
library(umap)

# Visualising data using UMAP
gbox.umap <- umap(gbox_filtered)

# Cell clustering of G-box related genes
colours = rainbow(length(unique(clust)))

plot(gbox.umap$layout[,1], 
     gbox.umap$layout[,2],
     col = colours[clust[includeCells]], 
     pch = 20,
     main = 'UMAP Flower', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2')


# PCA prior to UMAP:

pca = prcomp(gbox_filtered, 
             scale. = T, 
             rank. = 5)
gbox.pca.umap <- umap(pca$x)

colours = rainbow(length(unique(clust)))

plot(gbox.pca.umap$layout[,1], 
     gbox.pca.umap$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20, 
     main = 'PCA UMAP Flower', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2')






# Optional Task: 2024 Week 3 Homework - Gene Expression Within Clusters


# Previously scripted requirements for this task follow:

# install.packages('Matrix')
# library(Matrix)

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

# alternative filename='dev/utilities/dataprocessingHelperFunctions.R'
# source('dev/utilities/dataprocessingHelperFunctions.R')

# a=load('data/GSE226097_flower_230221.RData')
# araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)



# Average Expression of Each Gene in Each Cluster

clustAsNumbers = as.numeric(paste(clust))

geneExpByCluster = apply(gbox,
                         1,
                         function(i) {
                           sapply(0:(length(unique(clust))-1), 
                                  function(j) {
                                    ids = which(clustAsNumbers == j)
                                    mean(i[ids])
                                  })
                         })

colnames(geneExpByCluster) = rownames(gbox)

# Average gene expression of each gene across each cluster
dim(geneExpByCluster) 



# Visualising Large Tables as Heatmaps

install.packages('pheatmap') 
library('pheatmap')

# Plot gene expression changes as heatmap
pheatmap(geneExpByCluster, scale = 'column')

# Read in table of cell types and what clusters they belong to
clustLabs = read.table('data/clusterLabels.txt', header = T, sep = '\t')

# Meanings of clusters are unknown here, the clustLabs Organ col stores these
# Manually entered, so check data for cluster names
unique(clustLabs[,'Organ'])

# The organ analysed here is Flower. Below saves organ names and prints Flower clusters
organ='Flower'
simpleNames = clustLabs[which(clustLabs[,'Organ']==organ), "Cell.type.suggested"]
print(simpleNames)

# Set row names to those Flower cluster names
rownames(geneExpByCluster) = simpleNames

# Redo heatmap
pheatmap(geneExpByCluster, scale = 'column')

# Save underlying table as a file
write.table(geneExpByCluster, 'data/Flower_avgExpressionByCluster.txt')



# Focusing on Transcription Factors
# Same analysis, but focusing on the TFs that bind to perfect G-boxes

# A list of the Transcription Factors 
tfs = unique(araboxcis[,1])

# Filtered list for only TFs also included in the single-cell RNA-seq dataset
tfSubs = tfs[which(tfs %in% rownames(gbox))]

# Heatmap of the TFs that bind to perfect G-boxes
pheatmap(geneExpByCluster[, tfSubs], scale = 'column')



# What is the correlation between expression in TFs and potential downstream targets? 
# Does the TF regulate downstream gene?

# Correlations between every TF and every potential 
corMat = sapply(tfSubs, function(tf) {
  apply(geneExpByCluster, 
        2,
        function(gene){
          cor(geneExpByCluster[,tf], gene, method = 'spearman')
        })
})

# Dimensions of correlation testing
dim(corMat)

# Heatmap of the correlation testing
pheatmap(corMat)



# Looking at correlation between TF expression and downstream gene expression in individual cells

# TFs that are highly correlated with downstream genes when looking at cell type level
id = which(corMat > 0.8 & corMat != 1, arr.ind = TRUE)

# Dimensions of the above
dim(id)


# Neither Pearson's correlation or Spearman's are good choices, because there is a huge amount of 
#    zero-inflation when calculating correlation at an individual cell level. 
#    Instead, use the odds ratio of observing both the TF and downstream gene within the same cell.

# Odds ratio = 1 if no. of times two factors occur simultaneously is expected. 
# Odds ratio >1 if the two factors appear more 'positively' correlated with each other.
# Odds ratio <1 if the two factors are 'negatively' correlated with each other. 

# The following histogram plots the log-odds ratio. 
# Log of 1 = 0, log of <1 = negative, log >1 = positive. 

# Dimensions of correlation testing
dim(corMat)

# Dimensions of the TFs highly correlated with downstream genes
dim(id)

# Plotting the log-odds ratio
tVal = apply(id,
             1,
             function(i) {
               row = rownames(corMat)[i[1]]
               col = colnames(corMat)[i[2]]
               
               inBoth = length(which(gbox[row,] > 0 & gbox[col,] > 0))
               inNone = length(which(gbox[row,] == 0 & gbox[col,] == 0))
               inFirst = length(which(gbox[row,] > 0 & gbox[col,] == 0))
               inSecond = length(which(gbox[row,] == 0 & gbox[col,] > 0))
               
               matTemp = matrix(c(inBoth, inFirst, inSecond, inNone), ncol = 2)
               c(inBoth, inFirst, inSecond, inNone, (inBoth * inNone) / (inFirst * inSecond))
             })

hist(log(tVal[5,]), 
     main = 'Flower',
     xlab = 'Log Odds Ratio')



# Compare between datasets

thresh = exp(1)

idDoublePositive = which(tVal[5,] > thresh)
doublePositive = id[idDoublePositive,]
doublePositive[,1] = rownames(corMat)[doublePositive[,2]]
doublePositive[,2]  = colnames(corMat)[as.numeric(doublePositive[,2])]
doublePositive = cbind(doublePositive[,2], 
                       doublePositive[,1]) # TF comes before the target
doublePositive = cbind(doublePositive,
                       tVal[5,idDoublePositive])

# Save the file
write.table(doublePositive, 
            file = 'data/Flower_doublePositives.txt', 
            sep = '\t',
            row.names = F)

# Assembling a set of regulatory pairs that were positive in the first analysis, 
#    and negative in the second analysis. 

thresh = 1
idSimpson = which(tVal[5,] < thresh) 
Simpson = id[idSimpson,1]

Simpson[,1] = rownames(corMat)[Simpson[,1]]
Simpson[,2] = colnames(corMat)[as.numeric(Simpson[,2])]
Simpson = cbind(Simpson[,2], 
                Simpson[,1]) 
Simpson = cbind(Simpson,
                tVal[idSimpson])

# Save the file
write.table(Simpson,
            file = 'data/Flower_SimpsonPairs.txt',
            sep = '\t',
            row.names = F)





