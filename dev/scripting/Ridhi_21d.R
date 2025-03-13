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

##Week 3____________________________________________________________________________
library(GENIE3)
library(Matrix)
source('dev/utilities/dataprocessingHelperFunctions.R')

#Loading the dataset
a=load('data/GSE226097_rosette_21d_230221.RData')
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

#Making so that cluster is not be factors
clustAsNumbers=as.numeric(paste(clust))

geneExpByCluster = apply(gbox, 1, function(i){
  sapply(0:(length(unique(clust))-1), function(j){
    ids=which(clustAsNumbers==j)
    mean(i[ids])
  })
})

colnames(geneExpByCluster)=rownames(gbox)
dim(geneExpByCluster)

#Doing 5 tree network
net=GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees=5)

save(net, file='rosette_network_nTree_5.RData')

ginieOutput=convertToAdjacency(net, 0.05)
dim(ginieOutput)
ginieOutput[1:10,]

#Loading the network
a=load('rosette_network_nTree_5.RData')
newNet=GENIE3::getLinkList(net)
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

#Getting the set of unique genes in your new network
genesInNet=unique(c(newNet[,1], newNet[,2]))

#Filtering the AraBOXcis network to only contain genes that are in your new network
araboxcisFiltered=araboxcis[which(araboxcis[,1] %in% genesInNet & araboxcis[,2] %in% genesInNet),]

#Extracting the top edges in your new network, to make your network the same size as the araboxcisFiltered network.
newNetTopEdges=newNet[1:length(araboxcisFiltered[,1]),]

#Reformatting edges so it is more straightforward to compare them
edgesNew=paste(newNetTopEdges[,1], newNetTopEdges[,2], sep='_')
edgesOld=paste(araboxcisFiltered[,1], araboxcisFiltered[,2], sep='_')

#Coming up with all the different parts of a Venn Diagram:

#The overlap
length(which(edgesNew %in% edgesOld))

#In new network only
length(which(! (edgesNew %in% edgesOld)))

#In old network only
length(which(! (edgesOld %in% edgesNew)))

#Finding Important genes in Network
tfsNew=table(newNetTopEdges[,1])
tfsOld=table(araboxcisFiltered[,1])[names(tfsNew)]

#Histogram of Degrees
hist(as.numeric(tfsNew), main='SinceAraBOXcis', xlab='degree of TFs')

#Checking if the same TFs have high degrees in AraBOXcis and our new network:
plot(as.numeric(tfsNew), as.numeric(tfsOld), xlab='degree in SinceAraBOXcis', ylab='degree in AraBOXcis')

#Printing out the 20 TFs with highest degrees.  
sort(tfsNew, decreasing=TRUE)[1:20]

#Loading more libraries
library(igraph)
library(network)

#Going to heatmap
library(pheatmap)
simple_network <- graph_from_edgelist(as.matrix(newNetTopEdges[,c(1,2)]))

node_betweenness_all <- betweenness(simple_network)
node_betweenness=node_betweenness_all[which(node_betweenness_all>0)]
sort(node_betweenness, decreasing=TRUE)[1:20]

plot(sort(node_betweenness))

#abline(h=5000)
node_centrality_all <- alpha_centrality(simple_network)
node_centrality=node_centrality_all[which(node_centrality_all>0)]
sort(node_centrality, decreasing=TRUE)[1:20]

plot(sort(node_centrality))

#abline(h=5000)
node_hub_all <- hub_score(simple_network)$vector
node_hub=node_hub_all[which(node_hub_all>0)]
sort(node_hub, decreasing=TRUE)[1:20]

plot(sort(node_hub))

#abline(h=0.6)
plot(node_betweenness_all, node_centrality_all)
plot(node_hub_all, node_centrality_all)
plot(node_hub_all, node_betweenness_all)

#Associations between GO terms in the network
a=load('data/functionalData.RData')
source('dev/utilities/dataprocessingHelperFunctions.R')

pafwayOut=pafway(GOconcat, newNetTopEdges, unique(goOfInterest))
rownames(pafwayOut)=colnames(pafwayOut)

#Filtering to only include rows and columns with at least one significant factor:
atLeastOneSigRow=which(apply(pafwayOut, 1, function(i){length(which(i<0.05))})>0)

atLeastOneSigCol=which(apply(pafwayOut, 2, function(i){length(which(i<0.05))})>0)

pafwayInterestingOnly=pafwayOut[atLeastOneSigRow, atLeastOneSigCol]

#Here, the terms associated with the columns are upstream of the terms associated with the rows.  The values correspond to p-values, so smaller is more significant.
pheatmap(pafwayInterestingOnly)

#Re-doing zooming to the most significant associations by taking a log
pheatmap(log(pafwayInterestingOnly, 10))

##Doing the venn diagram task_________________________________________________
# Load required libraries
library(ggvenn)
library(VennDiagram)

# Defining  values
venn_data <- list(
  "New Network" = rep(1, 40072 + 9381),
  "Old Network" = rep(1, 40072 + 9381)
)


# Plotting the diagram 
venn.plot <- draw.pairwise.venn(
  area1 = 40072 + 9381,  
  area2 = 40072 + 9381,  
  cross.area = 9381,  
  category = c("New", "Old"),  # Shorter labels
  fill = c("lightpink", "lightblue"),
  alpha = 0.5,
  cat.cex = 1.2,  # Reduce category font size
  cex = 1.2      # Reduce numbers font size
)
