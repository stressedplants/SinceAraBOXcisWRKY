#install.packages('Matrix')
library(Matrix)

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#alternative filename='dev/utilities/dataprocessingHelperFunctions.R'
source('dev/utilities/dataprocessingHelperFunctions.R')

#Load the single cell data for flowers.
a=load('data/GSE226097_silique_230221.RData')  ########CHANGE THIS TO YOUR FILE!!!

#Load the original AraBOXcis network that was trained on bulk RNA-seq in seedlings 
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

# What is contained in this data set
print(a)

dim(gbox)

length(clust)

rownames(gbox)[1:10]

colnames(gbox)[1:10]

as.matrix(gbox[1:50, 1:3])

clust[1:8]

table(clust)

plot(table(sort(clust)), xlab='cluster name', ylab='number of cells', main='Rosette')

dim(araboxcis)

araboxcis[1:4,]

hist(araboxcis[,3])

tfs=unique(araboxcis[,1])

tfSubs=tfs[which(tfs %in% rownames(gbox))]

length(tfSubs)

# Filter genes and cells with very low values

dim(gbox)

thresh=0.01
numberGenesPerCell=apply(gbox, 2, function(i){length(which(i>0))})
includeCells=which(numberGenesPerCell>(thresh*dim(gbox)[1]))

gbox_filtered=gbox[,includeCells]

dim(gbox_filtered)

numberCellsPerGene=apply(gbox_filtered, 1, function(i){length(which(i>0))})
includeGenes=which(numberCellsPerGene>(thresh*dim(gbox_filtered)[2]))
gbox_filtered=gbox_filtered[includeGenes,]
dim(gbox_filtered)

# First visualisation with UMAP

install.packages('umap')
library(umap)

#Let us visualise it using UMAP
gbox.umap <- umap(gbox_filtered)

#Do the cell type clusters group together if we only look at G-box related genes?
colours=rainbow(length(unique(clust)))
plot(gbox.umap$layout[,1], gbox.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='UMAP Rosette', xlab='UMAP Component 1', ylab='UMAP Component 2')

pca=prcomp(gbox_filtered, scale.=T, rank.=5)
gbox.pca.umap <- umap(pca$x)

colours=rainbow(length(unique(clust)))
plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], col=colours[clust[includeCells]], pch=20, Rosette='PCA UMAP Rosette', xlab='UMAP Component 1', ylab='UMAP Component 2')



            








