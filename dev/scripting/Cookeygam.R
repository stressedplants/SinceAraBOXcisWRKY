# Week 1 Homework:  Downloading data (G-box) and doing some simple data visualisations.

# Installing packages(matrix) 
library(Matrix)

#if (!require("BiocManager", quitely = True))
# Install.packages("BioManager")
library(BiocManager)


# Loading datasets for seed_odays and AraBoxcis
# Loading single cell data sets ('Seed_0days')
Seed0Days = load("C:/Users/User/OneDrive/Documents/Seed_0days_Cookeygam/GSE226097_seed_0d_230221.RData")

# Loading the original AraBoxcis network that was trained on bulk RNA-seq in seedlings
AraBoxcis = read.csv("C:/Users/User/OneDrive/Documents/Seed_0days_Cookeygam/gboxNetwork22C.csv", header=TRUE)

# What is contained in this data set?
print(Seed0Days)
dim(gbox)
length(clust)
rownames(gbox) [1:10]
colnames(gbox) [1:10]
as.matrix(gbox[1:50, 1:3])
clust[1:8]
table(clust)

# Plotting graph using data from table
plot(table(sort(clust)), xlab = 'cluster name', ylab = 'number of cells', main = 'Seed0Days')

# Reading data from "Araboxcis"
dim(AraBoxcis)
AraBoxcis[1:4,]

# Making a Histogram of the scores to decide edges to be a part of the network or not
hist(AraBoxcis[,3])

# Getting a list of transcription factors
transcription_fs = unique(AraBoxcis[,1])

tfSubs = transcription_fs[which(transcription_fs %in% rownames(gbox))]
length(tfSubs)
dim(gbox)

# Getting rid of "Cells" that have less than 1% of the genes expressed
thresh = 0.01
NumberGenesPerCell = apply(gbox, 2, function(i){length(which(i>0))})
IncludeCells = which(NumberGenesPerCell > (thresh * dim(gbox)[1]))

gbox_filtered = gbox[,IncludeCells]
dim(gbox_filtered)

gbox[,IncludeCells]

# Getting rid of "Genes" that are expressed in less than 1% of the cells
NumberCellsPerGene = apply(gbox_filtered, 1, function(i){length(which(i>0))})
IncludeGenes = which(NumberCellsPerGene > (thresh * dim(gbox_filtered)[2]))

gbox_filtered = gbox_filtered[IncludeGenes,]
dim(gbox_filtered)

# Visualization with UMAP
install.packages('umap')
library(umap)

# Visualization using UMAP
gbox.umap = umap(gbox_filtered)

# G-box related gene
colours = rainbow(length(unique(clust)))
plot(gbox.umap$layout[,1], gbox.umap$layout[,2], col = colours[clust[IncludeCells]],
     pch = 20, main = 'UMAP Seed0Days', xlab = 'UMAP Component 1', ylab = 'UMAP Component 2')

# PCA UMAP Seed0Days'
PCA = prcomp(gbox_filtered, scale. = T, rank. = 5)
gbox.PCA.umap = umap(PCA$x)

# G-box related gene
colours = rainbow(length(unique(clust)))
plot(gbox.PCA.umap$layout[,1], gbox.PCA.umap$layout[,2], col = colours[clust[IncludeCells]],
     pch = 20, main = 'PCA UMAP Seed0Days', xlab = 'UMAP Component 1', ylab = 'UMAP Component 2')
