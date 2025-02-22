##Week 1 ________________________________________________________
#Loading the library
library(Matrix)

#Loading the helper function
source('dev/utilities/dataprocessingHelperFunctions.R')

#Loading the data
a=load('data/GSE226097_rosette_21d_230221.RData')

#Loading the original AraBOXcis network that was trained on bulk RNA-seq in seedlings 
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

##Contents of the dataset _______________________________________
#Viewing contents of 'a'
print(a)

#Checking dimentions of 'gbox'
dim(gbox)

#Checking the length of 'clust'
length(clust)

#Checking the  names of the first 10 rows
rownames(gbox)[1:10]

#Checking column names
colnames(gbox)[1:10]

#Viewing a 50X3 subset of the whole Matrix as a table
as.matrix(gbox[1:50, 1:3])

#Viewing the cluster designations of the first 8 cells
clust[1:8]

#Viewing the cell count in each cluster
table(clust)

#Plotting cell count against the clusters
plot(table(sort(clust)), xlab='cluster name', ylab='number of cells', main='Rosette')

#Checking dimentions of araboxcis matrix
dim(araboxcis)

#Viewing the first 4 rows of araboxcis
araboxcis[1:4,]

#Plotting a histogram of edge score frequencies
hist(araboxcis[,3])

#Filtering unique values from transcription factor column of araboxcis
tfs = unique(araboxcis[,1])

#Filtering tf, for transcription factors in scRNA data dataset
tfSubs=tfs[which(tfs %in% rownames(gbox))]

#Checking the lenght of unique TFs in 
length(tfSubs)

##Filtering cells and genes on the basis of very low values _________________
#Checking the dimentions of 'gbox'
dim(gbox)

#Filtering out cells with <1% genes expressed
thresh=0.01
numberGenesPerCell=apply(gbox, 2, function(i){length(which(i>0))})
includeCells=which(numberGenesPerCell>(thresh*dim(gbox)[1]))

gbox_filtered=gbox[,includeCells]

dim(gbox_filtered)

#Checking the dimensions of filtered gbox
numberCellsPerGene=apply(gbox_filtered, 1, function(i){length(which(i>0))})
includeGenes=which(numberCellsPerGene>(thresh*dim(gbox_filtered)[2]))
gbox_filtered=gbox_filtered[includeGenes,]
dim(gbox_filtered)

##First visualisation with UMAP ___________________________________________________
library(umap)

#Visualising the data using UMAP
gbox.umap <- umap(gbox_filtered)

#Visualising it using UMAP
gbox.umap <- umap(gbox_filtered)

#Checking cell type clustering under G-box related genes
colours=rainbow(length(unique(clust)))
plot(gbox.umap$layout[,1], gbox.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='UMAP Rosette', xlab='UMAP Component 1', ylab='UMAP Component 2')

#Plotting UMAP for PCA
pca=prcomp(gbox_filtered, scale.=T, rank.=5)
gbox.pca.umap <- umap(pca$x)

colours=rainbow(length(unique(clust)))
plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='PCA UMAP Rosette', xlab='UMAP Component 1', ylab='UMAP Component 2')
