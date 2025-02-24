# WEEK 1 69M PROJECT

# Loading library
library(Matrix)

# Loading helper function
source('HELPERFUNC/dataprocessingHelperFunctions.R')

# Loading the original AraBOXcis network trained on bulk RNA-seq in seedlings
araboxcis <- read.csv(file = 'data/gboxNetwork22C.csv',header = TRUE)

# Load the data 
seedling_6d= load("/Users/sakaoglu/Desktop/GSE226097_seedling_6d_230221.RData")


# The information of the data set-------------------------------------

# View contents of the data set 
print(seedling_6d)

# Check dimentions of 'gbox'
dim(gbox)

# Check the length of 'clust'
length(clust)

# Check the  names of the first 10 rows
rownames(gbox)[1:10]

# Check column names
colnames(gbox)[1:10]

# View a 50X3 subset of the whole Matrix as a table
as.matrix(gbox[1:50, 1:3])

# View the cluster designations of the first 8 cells
clust[1:8]

# View the cell count in each cluster
table(clust)

## Plot cell count against the clusters
plot(table(sort(clust)), xlab='cluster name', ylab='number of cells', main='GSE226097_seedling_6d_230221')

# Check dimentions of araboxcis matrix
dim(araboxcis)

# Ciew the first 4 rows of araboxcis
araboxcis[1:4,]

# Plot a histogram of edge score frequencies
hist(araboxcis[,3])

# Filter unique values from transcription factor column of araboxcis
tfs = unique(araboxcis[,1])

# Filter tf, for transcription factors in scRNA data dataset
tfSubs=tfs[which(tfs %in% rownames(gbox))]

# Check the lenght of unique TFs in 
length(tfSubs)

## Filter genes and cells -----------------------------

# Check the dimentions of 'gbox'
dim(gbox)

# Filter out cells with <1% genes expressed
thresh=0.01  

numberGenesPerCell = apply(gbox, 2, function(i){length(which(i>0))}) 
includeCells = which(numberGenesPerCell > (thresh*dim(gbox)[1])) 
gbox_filtered = gbox[,includeCells] 

# Check the dimensions of filtered gbox
dim(gbox_filtered)

# Filter out genes expressed in 1% of the cells
numberCellPerGene = apply(gbox, 1, function(i){length(which(i>0))})
includeGenes = which(numberCellPerGene > (thresh*dim(gbox)[2]))
gbox_filtered=gbox_filtered[includeGenes,]

dim(gbox_filtered)


## First visualisation with UMAP -------------------------------------------
library(umap)

# Visualising data using UMAP
gbox.umap <- umap(gbox_filtered)

# Checking cell clustering under G-box related genes
colours=rainbow(length(unique(clust)))

plot(gbox.umap$layout[,1], gbox.umap$layout[,2],
     col=colours[clust[includeCells]], pch=20,
     main='UMAP Seedling 6d', xlab='UMAP Component 1', ylab='UMAP Component 2')
dim(gbox_filtered)

# Checking cell clustering by PCA
pca=prcomp(gbox_filtered, scale.=T, rank.=5)
gbox.pca.umap <- umap(pca$x)

colours=rainbow(length(unique(clust)))
plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='PCA UMAP Seedling 6d', xlab='UMAP Component 1', ylab='UMAP Component 2')


