# WEEK 3 


# Loading libraries 
BiocManager::install("GENIE3")
library(GENIE3)

# Loading the data
seedling_6d= load("/Users/sakaoglu/Desktop/GSE226097_seedling_6d_230221.RData")

# Integrating network algorthim first try with 5 and move up
net = GENIE3(as.matrix(gbox),regulators = tfSubs, nTrees = 5) 
save(net, file='data/seedling-d6_network_nTree_5.RData') 

# Try to move up 
net8 = GENIE3(as.matrix(gbox),regulators = tfSubs, nTrees = 8) 
save(net8, file='data/seedling-d6_network_nTree_8.RData') 

# Try to move up
net10 = GENIE3(as.matrix(gbox),regulators = tfSubs, nTrees = 10) 
save(net10, file='data/seedling-d6_network_nTree_10.RData') 


#import the data in
load('data/seedling-d6_network_nTree_10.RData')
newNet = GENIE3::getLinkList(net10)
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

# Load the other data-----------------------------------------------------------
load('data/seedling-d6_network_nTree_5.RData')
load('data/seedling-d6_network_nTree_8.RData')
#-------------------------------------------------------------------------------

# Converting this to a table 
ginieOutput=convertToAdjacency(net10, 0.05)

# Looking at dimensions
dim(ginieOutput)

# The first 10 rows
ginieOutput[1:10,]
# A LOT OF NA VALUES WHY??


# Daphne's CODE
#gets the set of unique genes in your new network
genesInNet=unique(c(newNet[,1], newNet[,2]))

#filter the AraBOXcis network to only contain genes that are in your new network
araboxcisFiltered=araboxcis[which(araboxcis[,1] %in% genesInNet & araboxcis[,2] %in% genesInNet),]

#extract the top edges in your new network, to make your network the same size as the araboxcisFiltered network.
newNetTopEdges=newNet[1:length(araboxcisFiltered[,1]),]

#reformat edges so it is more straightforward to compare them
edgesNew=paste(newNetTopEdges[,1], newNetTopEdges[,2], sep='_')
edgesOld=paste(araboxcisFiltered[,1], araboxcisFiltered[,2], sep='_')

# Different parts of a Venn Diagram:

#The overlap
length(which(edgesNew %in% edgesOld))

#In NEW network only
length(which(! (edgesNew %in% edgesOld)))

#In OLD network only
length(which(! (edgesOld %in% edgesNew)))

#DO THE VENN DIAGRAM HERE-------------------------------------------------------

# Install packages 
install.packages("ggvenn")
install.packages("VennDiagram")

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
  fill = c("green", "lightblue"),
  alpha = 0.5,
  cat.cex = 1.2,  # Reduce category font size
  cex = 1.2,      # Reduce numbers font size
  fontface = "bold"
)
#-------------------------------------------------------------------------------

# Find important genes in the network 
tfsNew=table(newNetTopEdges[,1])
tfsOld=table(araboxcisFiltered[,1])[names(tfsNew)]

# Histogram of Tfs 
hist(as.numeric(tfsNew), main='SinceAraBOXcis', xlab='degree of TFs')

# Let's see if the same TFs have high degrees in AraBOXcis and our new network:
plot(as.numeric(tfsNew), as.numeric(tfsOld), xlab='degree in SinceAraBOXcis', ylab='degree in AraBOXcis')

#Let's print out the 20 TFs with highest degrees.  
sort(tfsNew, decreasing=TRUE)[1:20]

# EFP browser look at it--------------------------------------------------------




#-------------------------------------------------------------------------------

# Load packages 
install.packages("igraph")
library(igraph)
install.packages("network")
library(network)
library(pheatmap)

# Doing the network
simple_network <- graph_from_edgelist(as.matrix(newNetTopEdges[,c(1,2)]))

node_betweenness_all <- betweenness(simple_network)
node_betweenness=node_betweenness_all[which(node_betweenness_all>0)]
sort(node_betweenness, decreasing=TRUE)[1:20]

# Plot
plot(sort(node_betweenness))

# abline(h=5000)

node_centrality_all <- alpha_centrality(simple_network)
node_centrality=node_centrality_all[which(node_centrality_all>0)]
sort(node_centrality, decreasing=TRUE)[1:20]
# Out of memory ERROR??

# How to fix the error?
vcount(simple_network)  # Number of nodes
ecount(simple_network)  # Number of edges

# Other ways to fix the code:
node_centrality_all <- alpha_centrality(simple_network, alpha=0.9)

node_centrality=node_centrality_all[which(node_centrality_all>0)]
sort(node_centrality, decreasing=TRUE)[1:20]

plot(sort(node_centrality))

#abline(h=5000)

node_hub_all <- hub_score(simple_network)$vector
node_hub=node_hub_all[which(node_hub_all>0)]
sort(node_hub, decreasing=TRUE)[1:20]

plot(sort(node_hub))

#abline(h=0.6)

# Plot betweeness vs node centrality
plot(node_betweenness_all, node_centrality_all)

# Plot Hub all vs node centrality 
plot(node_hub_all, node_centrality_all)

# Plot Hub all vs betweeness 
plot(node_hub_all, node_betweenness_all)

# Which one is important? Understand whats going on 
# and discuss what is the most important graph from this?


#---------------------------------------------------------

# Find associations/ Occur in downstream or upstream?

# Upload the data
a=load('data/functionalData.RData')

pafwayOut=pafway(GOconcat, newNetTopEdges, unique(goOfInterest))
rownames(pafwayOut)=colnames(pafwayOut)

#Filter to only include rows and columns with at least one significant factor:
atLeastOneSigRow=which(apply(pafwayOut, 1, function(i){length(which(i<0.05))})>0)


atLeastOneSigCol=which(apply(pafwayOut, 2, function(i){length(which(i<0.05))})>0)

pafwayInterestingOnly=pafwayOut[atLeastOneSigRow, atLeastOneSigCol]

#Here, the terms associated with the columns are upstream of the terms associated with the rows.  The values correspond to p-values, so smaller is more significant.
pheatmap(pafwayInterestingOnly)

# Upgrade the image quality 
pheatmap(pafwayInterestingOnly, fontsize_row = 7, fontsize_col = 7)

#Let's re-do zooming to the most significant associations by taking a log
pheatmap(log(pafwayInterestingOnly, 10))

# Adjusting the image 
pheatmap(
  log10(pafwayInterestingOnly),  
  angle_col = 45,                
  fontsize_row = 6,              
  fontsize_col = 6,              
  width = 10,                    
  height = 10                      
)


#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
# WEEK 3 FOR nTREES= 5


#import the data in
load('data/seedling-d6_network_nTree_5.RData')
newNet = GENIE3::getLinkList(net)
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

#-------------------------------------------------------------------------------

# Converting this to a table 
ginieOutput=convertToAdjacency(net, 0.05)

# Looking at dimensions
dim(ginieOutput)

# The first 10 rows
ginieOutput[1:10,]
# A LOT OF NA VALUES WHY??


# Daphne's CODE
#gets the set of unique genes in your new network
genesInNet=unique(c(newNet[,1], newNet[,2]))

#filter the AraBOXcis network to only contain genes that are in your new network
araboxcisFiltered=araboxcis[which(araboxcis[,1] %in% genesInNet & araboxcis[,2] %in% genesInNet),]

#extract the top edges in your new network, to make your network the same size as the araboxcisFiltered network.
newNetTopEdges=newNet[1:length(araboxcisFiltered[,1]),]

#reformat edges so it is more straightforward to compare them
edgesNew=paste(newNetTopEdges[,1], newNetTopEdges[,2], sep='_')
edgesOld=paste(araboxcisFiltered[,1], araboxcisFiltered[,2], sep='_')

# Different parts of a Venn Diagram:

#The overlap
length(which(edgesNew %in% edgesOld))

#In NEW network only
length(which(! (edgesNew %in% edgesOld)))

#In OLD network only
length(which(! (edgesOld %in% edgesNew)))

#DO THE VENN DIAGRAM HERE-------------------------------------------------------

# Install packages 
install.packages("ggvenn")
install.packages("VennDiagram")

# Load required libraries
library(ggvenn)
library(VennDiagram)

# Defining  values
venn_data <- list(
  "New Network" = rep(1, 40140 + 9313),
  "Old Network" = rep(1, 40140 + 9313)
)

# Plotting the diagram 
venn.plot <- draw.pairwise.venn(
  area1 = 40140 + 9313,  
  area2 = 40140 + 9313,  
  cross.area = 9313,  
  category = c("New", "Old"),  # Shorter labels
  fill = c("green", "lightblue"),
  alpha = 0.5,
  cat.cex = 1.2,  # Reduce category font size
  cex = 1.2,      # Reduce numbers font size
  fontface = "bold"
)
#-------------------------------------------------------------------------------

# Find important genes in the network 
tfsNew=table(newNetTopEdges[,1])
tfsOld=table(araboxcisFiltered[,1])[names(tfsNew)]

# Histogram of Tfs 
hist(as.numeric(tfsNew), main='SinceAraBOXcis', xlab='degree of TFs')

# Let's see if the same TFs have high degrees in AraBOXcis and our new network:
plot(as.numeric(tfsNew), as.numeric(tfsOld), xlab='degree in SinceAraBOXcis', ylab='degree in AraBOXcis')

#Let's print out the 20 TFs with highest degrees.  
sort(tfsNew, decreasing=TRUE)[1:20]

# EFP browser look at it--------------------------------------------------------



#-------------------------------------------------------------------------------

# Load packages 
install.packages("igraph")
library(igraph)
install.packages("network")
library(network)
library(pheatmap)

# Doing the network
simple_network <- graph_from_edgelist(as.matrix(newNetTopEdges[,c(1,2)]))

node_betweenness_all <- betweenness(simple_network)
node_betweenness=node_betweenness_all[which(node_betweenness_all>0)]
sort(node_betweenness, decreasing=TRUE)[1:20]

# Plot
plot(sort(node_betweenness))

# abline(h=5000)

node_centrality_all <- alpha_centrality(simple_network)
node_centrality=node_centrality_all[which(node_centrality_all>0)]
sort(node_centrality, decreasing=TRUE)[1:20]
# Out of memory ERROR??

# How to fix the error?
vcount(simple_network)  # Number of nodes
ecount(simple_network)  # Number of edges

# Other ways to fix the code:
node_centrality_all <- alpha_centrality(simple_network, alpha=0.9)

node_centrality=node_centrality_all[which(node_centrality_all>0)]
sort(node_centrality, decreasing=TRUE)[1:20]

plot(sort(node_centrality))

#abline(h=5000)

node_hub_all <- hub_score(simple_network)$vector
node_hub=node_hub_all[which(node_hub_all>0)]
sort(node_hub, decreasing=TRUE)[1:20]

plot(sort(node_hub))

#abline(h=0.6)

# Plot betweeness vs node centrality
plot(node_betweenness_all, node_centrality_all)

# Plot Hub all vs node centrality 
plot(node_hub_all, node_centrality_all)

# Plot Hub all vs betweeness 
plot(node_hub_all, node_betweenness_all)

# Which one is important? Understand whats going on 
# and discuss what is the most important graph from this


#---------------------------------------------------------

# Find associations/ Occur in downstream or upstream?

# Upload the data
a=load('data/functionalData.RData')

pafwayOut=pafway(GOconcat, newNetTopEdges, unique(goOfInterest))
rownames(pafwayOut)=colnames(pafwayOut)

#Filter to only include rows and columns with at least one significant factor:
atLeastOneSigRow=which(apply(pafwayOut, 1, function(i){length(which(i<0.05))})>0)


atLeastOneSigCol=which(apply(pafwayOut, 2, function(i){length(which(i<0.05))})>0)

pafwayInterestingOnly=pafwayOut[atLeastOneSigRow, atLeastOneSigCol]

#Here, the terms associated with the columns are upstream of the terms associated with the rows.  The values correspond to p-values, so smaller is more significant.
pheatmap(pafwayInterestingOnly)

# Upgrade the image quality 
pheatmap(pafwayInterestingOnly, fontsize_row = 7, fontsize_col = 7)

#Let's re-do zooming to the most significant associations by taking a log
pheatmap(log(pafwayInterestingOnly, 10))

# Adjusting the image 
pheatmap(
  log10(pafwayInterestingOnly),  
  angle_col = 45,                
  fontsize_row = 6,              
  fontsize_col = 6,              
  width = 10,                    
  height = 10                      
)
