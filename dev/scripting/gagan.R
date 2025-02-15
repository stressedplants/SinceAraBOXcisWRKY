
# Week 1 ------------------------------------------------------------------
## Load data ---------------------------------------------------------------

# Loading libraries
library(Matrix)

# Loading helper functions
source('dev/utilities/dataprocessingHelperFunctions.R')

# Loading the data
a=load('data/GSE226097_seedling_12d_230221.RData')

# Loading the original AraBOXcis network trained on bulk RNA-seq in seedlings
araboxcis <- read.csv(file = 'data/gboxNetwork22C.csv',header = TRUE)


## What is contained in this data set? -------------------------------------

# View contents of 'a'
print(a)

# check dimentions of 'gbox'
dim(gbox)

# check the length of 'clust'
length(clust)


# check the  names of the first 10 rows
rownames(gbox)[1:10]

# check column names
colnames(gbox)[1:10]

# view a 50X3 subset of the whole Matrix as a table
as.matrix(gbox[1:50, 1:3])

# view the cluster designations of the first 8 cells
clust[1:8]

# view the cell count in each cluster
table(clust)

## plot cell count against the clusters
plot(table(sort(clust)), xlab='cluster name', ylab='number of cells', main='GSE226097_seedling_12d_230221')

# check dimentions of araboxcis matrix
dim(araboxcis)

# view the first 4 rows of araboxcis
araboxcis[1:4,]

# plot a histogram of edge score frequencies
hist(araboxcis[,3])

# filter unique values from transcription factor column of araboxcis
tfs = unique(araboxcis[,1])

# filter tf, for transcription factors in scRNA data dataset
tfSubs=tfs[which(tfs %in% rownames(gbox))]
        
# check the lenght of unique TFs in 
length(tfSubs)


## Filter genes and cells with very low values -----------------------------

# check the dimentions of 'gbox'
dim(gbox)



# filter out cells with <1% genes expressed
thresh=0.01   #1 % threshold

numberGenesPerCell = apply(gbox, 2, function(i){length(which(i>0))}) # make a table of Cell-IDs and genes per cell
includeCells = which(numberGenesPerCell > (thresh*dim(gbox)[1])) # set the threshold value of 1% of length of genes in gbox, and filter numberGenesPerCell

gbox_filtered = gbox[,includeCells] # filter 'gbox' using the includeCells as a filter of the columns

# check the dimensions of filtered gbox
dim(gbox_filtered)


# filter out genes expressed in 1% of the cells
numberCellPerGene = apply(gbox, 1, function(i){length(which(i>0))})
includeGenes = which(numberCellPerGene > (thresh*dim(gbox)[2]))
gbox_filtered=gbox_filtered[includeGenes,]

dim(gbox_filtered)




## First visualisation with UMAP -------------------------------------------
library(umap)

# visualising data using UMAP
gbox.umap <- umap(gbox_filtered)

# checking cell clustering under G-box related genes
colours=rainbow(length(unique(clust)))

plot(gbox.umap$layout[,1], gbox.umap$layout[,2],
     col=colours[clust[includeCells]], pch=20,
     main='GSE226097_seedling_12d_230221', xlab='UMAP Component 1', ylab='UMAP Component 2')

# UMAP of PCA
pca = prcomp(gbox_filtered, scale. = T, rank. = 5)
gbox.pca.umap <- umap(pca$x)

colours=rainbow(length(unique(clust)))
plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], 
     col=colours[clust[includeCells]], pch=20, 
     main='GSE226097_seedling_12d_230221_PCA', xlab='UMAP Component 1', ylab='UMAP Component 2')
