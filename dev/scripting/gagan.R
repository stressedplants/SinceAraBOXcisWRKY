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




# Week 1 Optional ---------------------------------------------------------

## Load packages and data --------------------------------------------------
library(Matrix)
source('dev/utilities/dataprocessingHelperFunctions.R')
a=load('data/GSE226097_seedling_12d_230221.RData')
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

## Average expression of each gene in each cluster -------------------------
#make clust not be factors
clustAsNumbers=as.numeric(paste(clust))

geneExpByCluster = apply(gbox, 1, function(i){
  sapply(0:(length(unique(clust))-1), function(j){
    ids=which(clustAsNumbers==j)
    mean(i[ids])
  })
})

colnames(geneExpByCluster)=rownames(gbox)
dim(geneExpByCluster)


## Visualising large tables as heatmaps ------------------------------------
library('pheatmap')

#pheatmap(geneExpByCluster, scale='column')

clustLabs=read.table('data/clusterLabels.txt', header=T, sep='\t') #read table of cell type in clusters

unique(clustLabs[,'Organ']) # view the names of the dev stage(organs) 

organ='Seedlings_12d' #set the 'organ' name for my data set: Seedlings_12d
simpleNames=clustLabs[which(clustLabs[,'Organ']==organ), "Cell.type.suggested"] #filter the cell types for the selected stage
print(simpleNames)

rownames(geneExpByCluster)=simpleNames # set the row names from the filtered list

pheatmap(geneExpByCluster, scale='column') # redo heatmap

write.table(geneExpByCluster, 'data/Seedlingd12_avgExpressionByCluster.txt') # write the table of average expression per cluster

## Focussing on transcription factors --------------------------------------
tfs=unique(araboxcis[,1]) #Select the unique elements 'from' column with TFs
tfSubs=tfs[which(tfs %in% rownames(gbox))] #select the tfs present in the gbox dataset
pheatmap(geneExpByCluster[,tfSubs], scale='column') 


## correlation b/w TF expression and  downstream targets -------------------
corMat=sapply(tfSubs, function(tf){
  apply(geneExpByCluster, 2, function(gene){
    cor(geneExpByCluster[,tf], gene, method='spearman')
  })
})

dim(corMat)

pheatmap(corMat)


## correlation b/w TF and targets per cell ---------------------------------
id=which(corMat>0.8 & corMat!=1, arr.ind = TRUE)
dim(id)

dim(corMat)

tVal=apply(id, 1, function(i){
  row=rownames(corMat)[i[1]]
  col=colnames(corMat)[i[2]]
  
  inBoth=length(which(gbox[row,]>0 & gbox[col,]>0))
  inNone=length(which(gbox[row,]==0 & gbox[col,]==0))
  inFirst=length(which(gbox[row,]>0 & gbox[col,]==0))
  inSecond=length(which(gbox[row,]==0 & gbox[col,]>0))
  matTemp=matrix(c(inBoth, inFirst, inSecond, inNone), ncol=2)
  c(inBoth, inFirst, inSecond, inNone, (inBoth*inNone)/(inFirst*inSecond))
})

hist(log(tVal[5,]), main='Flower', xlab='log odds ratio')

## Comparing between the data sets -----------------------------------------

#lists of genes with positive on cell type level & single cell level
thresh=exp(1)
idDoublePositive=which(tVal[5,]>thresh)
doublePositive=id[idDoublePositive,]
doublePositive[,1]=rownames(corMat)[doublePositive[,1]]
doublePositive[,2]=colnames(corMat)[as.numeric(doublePositive[,2])]
doublePositive=cbind(doublePositive[,2], doublePositive[,1]) #so TF comes before target
doublePositive=cbind(doublePositive, tVal[5,idDoublePositive])

#save file
write.table(doublePositive, file='data/seedling-d12_doublePositives.txt', sep='\t', row.names=F)

#list of regulatory pairs that were positive in the first analysis and negative in the second analysis
thresh=1
idSimpson=which(tVal[5,]<thresh)
Simpson=id[idSimpson,]

Simpson[,1]=rownames(corMat)[Simpson[,1]]
Simpson[,2]=colnames(corMat)[as.numeric(Simpson[,2])]
Simpson=cbind(Simpson[,2], Simpson[,1]) #so TF comes before target
Simpson=cbind(Simpson, tVal[idSimpson])

#save file
write.table(Simpson, file='data/seedling-d12_SimpsonPairs.txt', sep='\t', row.names=F)