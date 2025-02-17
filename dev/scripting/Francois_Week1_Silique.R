#Silique-week 1

#install.packages('Matrix')
#load libraries
library (tidyverse)
library(Matrix)
library(ggplot2)

#install packages
install.packages("BiocManager")

#install HelperFunctions.R
source('dev/utilities/dataprocessingHelperFunctions.R')

#load Silique data
a=load('data/GSE226097_silique_230221.RData')

#Load the original AraBOXcis network  trained on bulk RNA-seq in seedlings 
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

#viewing contents of data set
print(a)

#view dimension of gbox
dim(gbox)

#view no. of elements per cluster 
length(clust)

#view first 10 rownames of gbox
rownames(gbox)[1:10]

#view first 10 column names of box
colnames(gbox)[1:10]

#Extracting sub-portion(first 50 rows from first 3 columns)
as.matrix(gbox[1:50, 1:3])

#print cluster designations of first 8 cells
clust[1:8]

#use table to count no of cells per cluster
table(clust)

#plot frequency table 
plot(table(sort(clust)), xlab='cluster name', ylab='number of cells', main='Silique')

#checking data for gboxNetwork22C
dim(araboxcis)

#view first 4 rows of gbox
araboxcis[1:4,]

#plot histogram for araboxcis
hist(araboxcis[,3])

#View no of transcription factors
tfs=unique(araboxcis[,1])

#filter for transcription factors also present in scRNA seq data
tfSubs=tfs[which(tfs %in% rownames(gbox))]

#checking for transcription factors also present in scRNA-seq data
length(tfSubs)

#Check dimensions of gbox
dim(gbox)

#Filtering for cells with <1% genes expressed
thresh=0.01
numberGenesPerCell=apply(gbox, 2, function(i){length(which(i>0))})
includeCells=which(numberGenesPerCell>(thresh*dim(gbox)[1]))

gbox_filtered=gbox[,includeCells]

dim(gbox_filtered)

#Remove cells with <1% genes expressed
numberCellsPerGene=apply(gbox_filtered, 1, function(i){length(which(i>0))})
includeGenes=which(numberCellsPerGene>(thresh*dim(gbox_filtered)[2]))
gbox_filtered=gbox_filtered[includeGenes,]
dim(gbox_filtered)

#visualisation- install.packages('umap')
library(umap)

#visualisation using umap
gbox.umap <- umap(gbox_filtered)

#Do the cell type clusters group together if we only look at G-box related genes?
colours=rainbow(length(unique(clust)))
plot(gbox.umap$layout[,1], gbox.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='UMAP Silique', xlab='UMAP Component 1', ylab='UMAP Component 2')

#Plotting PCA
pca=prcomp(gbox_filtered, scale.=T, rank.=5)
gbox.pca.umap <- umap(pca$x)

#Plotting PCA UMAP
colours=rainbow(length(unique(clust)))
plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='PCA UMAP Silique', xlab='UMAP Component 1', ylab='UMAP Component 2')
