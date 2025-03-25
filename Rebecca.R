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


###----------------------------_________________________________________________

# WEEK 3

library(GENIE3)

net=GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees=5) #use more than 5 trees if your computer can handle it

save(net, file='rosette_network_nTree_5.RData')

# attempting 10 trees
net=GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees=5) #use more than 5 trees if your computer can handle it

save(net, file='rosette_network_nTree_5.RData')

# To convert data into table

load("rosette_network_nTree_5.RData")

ginieOutput=convertToAdjacency(net, 0.05)

### Look at the dimensions and look at first 10 edges in the network

dim(ginieOutput)

ginieOutput[1:10,]

### Load the network

a=load('rosette_network_nTree_5.RData') ### Remember to actually load in your specific network 
newNet=GENIE3::getLinkList(net)
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

### Find overlaps between the single cell network and the AraBOXcis network

#gets the set of unique genes in your new network
genesInNet=unique(c(newNet[,1], newNet[,2]))

#filter the AraBOXcis network to only contain genes that are in your new network
araboxcisFiltered=araboxcis[which(araboxcis[,1] %in% genesInNet & araboxcis[,2] %in% genesInNet),]

#extract the top edges in your new network, to make your network the same size as the araboxcisFiltered network.
newNetTopEdges=newNet[1:length(araboxcisFiltered[,1]),]

#reformat edges so it is more straightforward to compare them
edgesNew=paste(newNetTopEdges[,1], newNetTopEdges[,2], sep='_')
edgesOld=paste(araboxcisFiltered[,1], araboxcisFiltered[,2], sep='_')

#Now, you can come up with all the different parts of a Venn Diagram:

#The overlap
length(which(edgesNew %in% edgesOld))

#In new network only
length(which(! (edgesNew %in% edgesOld)))

#In old network only
length(which(! (edgesOld %in% edgesNew)))

### Venn diagram 
#DO THE VENN DIAGRAM HERE-----------------------

# Install packages 
install.packages("ggvenn")
install.packages("VennDiagram")

# Load required libraries
library(ggvenn)
library(VennDiagram)

# Defining  values
venn_data <- list(
  "New Network" = rep(1, 37787 + 8947),
  "Old Network" = rep(1, 37787 + 8947)
)

# Plotting the diagram 
venn.plot <- draw.pairwise.venn(
  area1 = 37787 + 8947,  
  area2 = 37787 + 8947,  
  cross.area = 8947,  
  category = c("New", "Old"),  # Shorter labels
  fill = c("green", "lightblue"),
  alpha = 0.5,
  cat.cex = 1.2,  # Reduce category font size
  cex = 1.2,      # Reduce numbers font size
  fontface = "bold"
)

### Find important genes in network

tfsNew=table(newNetTopEdges[,1])
tfsOld=table(araboxcisFiltered[,1])[names(tfsNew)]

#histogram of degrees should look like an exponential distribution, because biological networks are a kind of network called 'scale-free', meaning that most TFs only regulate a small number of genes, but a few are influential hubs.

hist(as.numeric(tfsNew), main='SinceAraBOXcis', xlab='degree of TFs')

#Let's see if the same TFs have high degrees in AraBOXcis and our new network:
plot(as.numeric(tfsNew), as.numeric(tfsOld), xlab='degree in SinceAraBOXcis', ylab='degree in AraBOXcis')

#Let's print out the 20 TFs with highest degrees.  
sort(tfsNew, decreasing=TRUE)[1:20]

### Calculate all the different metrics of TF importance

install.packages("igraph")
library(igraph)

install.packages("network")
library(network)

library(pheatmap)
simple_network <- graph_from_edgelist(as.matrix(newNetTopEdges[,c(1,2)]))

node_betweenness_all <- betweenness(simple_network)
node_betweenness=node_betweenness_all[which(node_betweenness_all>0)]
sort(node_betweenness, decreasing=TRUE)[1:10]

plot(sort(node_betweenness))

#abline(h=5000)

node_centrality_all <- alpha_centrality(simple_network, alpha=0.9)

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

############
############ Find associations between GO terms in the network

a=load('data/functionalData.RData')
source('dev/utilities/dataprocessingHelperFunctions.R')

pafwayOut=pafway(GOconcat, newNetTopEdges, unique(goOfInterest))
rownames(pafwayOut)=colnames(pafwayOut)

#Filter to only include rows and columns with at least one significant factor:
atLeastOneSigRow=which(apply(pafwayOut, 1, function(i){length(which(i<0.05))})>0)


atLeastOneSigCol=which(apply(pafwayOut, 2, function(i){length(which(i<0.05))})>0)

pafwayInterestingOnly=pafwayOut[atLeastOneSigRow, atLeastOneSigCol]

#Here, the terms associated with the columns are upstream of the terms associated with the rows.  The values correspond to p-values, so smaller is more significant.
pheatmap(pafwayInterestingOnly, fontsize_row = 8, fontsize_col = 8)

#Let's re-do zooming to the most significant associations by taking a log
pheatmap(
  log10(pafwayInterestingOnly),  
  angle_col = 45,                
  fontsize_row = 6,              
  fontsize_col = 6,              
  width = 10,                    
  height = 10)



save(node_betweenness, file = 'data/centrality_rosette30.RData')


save(node_betweenness_all, 
     file = 'data/centrality_all_rosette30.RData')

save(node_betweenness_all,node_centrality_all,node_hub_all, 
     file = 'data/centrality_all_rosette30.RData')



#########################################################################################

# Open a new file in the directory named exported_data

# 1. Export your network edge list
# This uses your GENIE3 network 
write.csv(newNetTopEdges, file = "exported_data/Rosette_30d_network.csv", row.names = FALSE)

# 2. Export your GO enrichment results (the pafway output)
# First create a long-format table from your pafwayInterestingOnly matrix
go_results <- data.frame()
for (row_term in rownames(pafwayInterestingOnly)) {
  for (col_term in colnames(pafwayInterestingOnly)) {
    go_results <- rbind(go_results, data.frame(
      DownstreamTerm = row_term,
      UpstreamTerm = col_term,
      pValue = pafwayInterestingOnly[row_term, col_term]
    ))
  }
}

# Export the GO results
write.csv(go_results, file = "exported_data/Rosette_30d_GO_results.csv", row.names = FALSE)

# 3. Export the GO terms that were found to be significant 
sig_go_terms <- unique(c(rownames(pafwayInterestingOnly), colnames(pafwayInterestingOnly)))
write.csv(data.frame(GOTerm = sig_go_terms), file = "exported_data/Rosette_30d_sig_GO_terms.csv", row.names = FALSE)

# Optional: Export the full pafwayOut matrix before filtering for significant terms
# This contains all GO terms analyzed, not just the significant ones
write.csv(pafwayOut, file = "exported_data/Seedling_6d_full_pafway_matrix.csv")



#########################______________________________________________________


# Creating the network in cytoscape
# Select top 1000 edges
top_edges <- newNetTopEdges[1:1000,]

write.table(top_edges, file="exported_data/rosette_30d_top1000_edges.csv", sep=",", row.names=FALSE, quote=FALSE)
