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



# Week 3 ------------------------------------------------------------------


### If reload required ------------------------------------------------------
# # Loading helper functions
# source('dev/utilities/dataprocessingHelperFunctions.R')
# 
# # Loading the data
# a=load('data/GSE226097_seedling_12d_230221.RData')
# 
# # Loading the original AraBOXcis network trained on bulk RNA-seq in seedlings
# araboxcis <- read.csv(file = 'data/gboxNetwork22C.csv',header = TRUE)
# 
# # filter unique values from transcription factor column of araboxcis
# tfs = unique(araboxcis[,1])
# 
# # filter tf, for transcription factors in scRNA data dataset
# tfSubs=tfs[which(tfs %in% rownames(gbox))]



## Creating the network ----------------------------------------------------

library(GENIE3) #load the GENIE3 to create the network with

library(tictoc) #load the tictoc library to measure the runtime 
library(beepr) #load the beepr library to beep when code has been run

tic('nTree = 11, nCore = 4') #starting the timer

net = GENIE3(as.matrix(gbox),regulators = tfSubs, nTrees = 11, nCores = 4) #creating the network (network sizes=5)
save(net, file='data/seedling-d12_network_nTree_11_nCore_4.RData') #save the network

### remember to change saved file name!!!!!!!!!!
toc() #timer end
beep(3) #play sound at the end of the code

# convert the network output into an adjecency matrix
ginieOutput=convertToAdjacency(net, 0.05)
dim(ginieOutput)
ginieOutput[1:10,]

# loading the network
load('data/seedling-d12_network_nTree_11_nCore_4.RData')
newNet = GENIE3::getLinkList(net)
araboxcis = read.csv('data/gboxNetwork22C.csv', header = T)


## Overlaps b/w networks ---------------------------------------------------

# create a set of unique genes in the new network
genesInNet = unique(c(newNet[,1],newNet[,2]))

# filter the araboxcis dataset to only include genes from the new network
araboxcisFiltered = araboxcis[which(araboxcis[,1] %in% genesInNet & araboxcis[,2] %in% genesInNet),]

# extract the top edges of the network,to make it the same size as araboxcisFiltered network
newNetTopEdges=newNet[1:length(araboxcisFiltered[,1]),]

# reformat edges for easier comparision
edgesNew=paste(newNetTopEdges[,1], newNetTopEdges[,2], sep='_')
edgesOld=paste(araboxcisFiltered[,1], araboxcisFiltered[,2], sep='_')

# create a Venn diagram of overlap
library(ggplot2)
library(ggvenn)

d = data_frame(edgesOld,edgesNew)

all = unique(c(d$edgesNew,d$edgesOld))

ven = tibble(edge = all,
             old_network = all %in%  d$edgesOld,
             new_network = all %in%  d$edgesNew)

ggplot(ven, aes(A = old_network, B = new_network)) + 
  geom_venn(fill_color = c('green','cyan')) +
  coord_fixed() + theme_void()


## important genes in network ----------------------------------------------

tfsNew = table(newNetTopEdges[,1]) # list the tfs from the new network and their frequency of occurrence

tfsOld = table(araboxcisFiltered[,1])[names(tfsNew)] # list the tfs from the new network and their frequency of occurrence and sort it according to tfsNew

hist(as.numeric(tfsNew), main='SinceAraBOXcis', xlab='degree of TFs') # view a friquency distribution of degree of TFs

# compare the degree of TFs in new and old network
plot(as.numeric(tfsNew), as.numeric(tfsOld), xlab='degree in SinceAraBOXcis', ylab='degree in AraBOXcis')

# print the top 20 TFs in new network with the highest degrees 
sort(tfsNew, decreasing=TRUE)[1:20]

# calculate different metricies of TF importance
library(igraph)
library(network)

library(pheatmap)
## node betweenness
simple_network <- graph_from_edgelist(as.matrix(newNetTopEdges[,c(1,2)]))

node_betweenness_all <- betweenness(simple_network)
node_betweenness <- node_betweenness_all[which(node_betweenness_all>0)]
sort(node_betweenness, decreasing = T)[1:20]

plot(sort(node_betweenness))

## node centrality !!!! DIDN'T WORK
#abline(h=5000)

node_centrality_all <- alpha_centrality(simple_network)
gcnode_centrality=node_centrality_all[which(node_centrality_all>0)]
sort(node_centrality, decreasing=TRUE)[1:20]

plot(sort(node_centrality))

## node hub
node_hub_all <- hub_score(simple_network)$vector
node_hub=node_hub_all[which(node_hub_all>0)]
sort(node_hub, decreasing=TRUE)[1:20]

plot(sort(node_hub))

# plot betweenness against centrality

plot(node_betweenness_all, node_centrality_all)

# plot hub against centrality

plot(node_hub_all, node_centrality_all)

# plot hubs against betweenness

plot(node_hub_all, node_betweenness_all)






## GO analysis of the network  ---------------------------------------------

a = load('data/functionalData.RData')
source('dev/utilities/dataprocessingHelperFunctions.R')

pafwayOut=pafway(GOconcat, newNetTopEdges, unique(goOfInterest))
rownames(pafwayOut)=colnames(pafwayOut)

#Filter to only include rows and columns with at least one significant factor:
atLeastOneSigRow=which(apply(pafwayOut, 1, function(i){length(which(i<0.05))})>0)

atLeastOneSigCol=which(apply(pafwayOut, 2, function(i){length(which(i<0.05))})>0)

pafwayInterestingOnly=pafwayOut[atLeastOneSigRow, atLeastOneSigCol]

#Here, the terms associated with the columns are upstream of the terms associated with the rows.
#The values correspond to p-values, so smaller is more significant.
pheatmap(pafwayInterestingOnly, main = 'Heatmap of Significancy, Gene Ontology',
         color=colorRampPalette(c("hotpink","lightyellow", "cyan"))(50))

#Let's re-do zooming to the most significant associations by taking a log
pheatmap(log(pafwayInterestingOnly, 10),main = 'Heatmap of Log Significance, Gene Ontology')









