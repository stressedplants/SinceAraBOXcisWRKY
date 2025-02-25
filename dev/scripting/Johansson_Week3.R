

# Week 3 Homework - Making and Exploring a Gene Network


# Necessary details from Week 1 homework

# Matrix package
library(Matrix)

# Loading helper functions
source('dev/utilities/dataprocessingHelperFunctions.R')

# if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

# Loading the flower data
a = load("~/Documents/Bioinformatics/Araboxcis/data/flowerdata.RData", verbose=TRUE)

# Loading the original AraBOXcis network trained on bulk RNA-seq in seedlings
araboxcis <- read.csv(file = 'data/gboxNetwork22C.csv',header = TRUE)

# A list of the Transcription Factors 
tfs = unique(araboxcis[,1])

# Filtered list for only TFs also included in the single-cell RNA-seq dataset
tfSubs = tfs[which(tfs %in% rownames(gbox))]



# Making a network from scratch

BiocManager::install("GENIE3")
library(GENIE3)

# GENIE3 = network inference algorithm. Its parameters include:
# 1. A gene expression dataset (gbox in this case)
# 2. A list of the regulators (tfSubs)
# 3. Various parameters related to how the random forest runs. 
#    We will decrease the default no. of trees in the random forest to speed up the code. 

# With 5 trees
net = GENIE3(as.matrix(gbox), 
             regulators = tfSubs, 
             nTrees=5)

# Save the tree(s!)
save(net,
    file = 'Flower_Network_nTree_5.RData')



# Converting the tree to a table
ginieOutput = convertToAdjacency(net, 0.05)

# Dimensions of the network
dim(ginieOutput)

# First 10 edges of the network
ginieOutput[1:10,]



# Load the network
a = load('Flower_Network_nTree_5.RData')
newNet = GENIE3::getLinkList(net)
araboxcis = read.csv('data/gboxNetwork22C.csv', header = T)



# Find overlaps between the single-cell gene network and the AraBOXcis network

# Gets the set of unique genes in the network
genesInNet = unique(c(newNet[,1], newNet[,2]))

# Filtered AraBOXcis network to only contain genes in the new network 
araboxcisFiltered = araboxcis[which(araboxcis[,1] %in% genesInNet & araboxcis[,2] %in% genesInNet),]

# Extracting top edges in the new network, making the network the same size as the araboxcisFiltered network
newNetTopEdges = newNet[1:length(araboxcisFiltered[,1]),]

# Reformatting edges to make comparisons easier
edgesNew = paste(newNetTopEdges[,1], newNetTopEdges[,2], sep = '_')
edgesOld = paste(araboxcisFiltered[,1], araboxcisFiltered[,2], sep = '_')

# Overlap
length(which(edgesNew %in% edgesOld))

# In the new network only
length(which(! (edgesNew %in% edgesOld)))

# In the old network only
length(which(! (edgesOld %in% edgesNew)))



# Find important genes in the network

# Different methods include finding:
# - Well-connected hubs (hub_score).
# - Influential/central genes (alpha centrality).
# - Genes that bridge between different regions in the network (betweenness). 
# The following calculates all three. 


# Degree = TFs that have the most edges.

# Transcription Factor edges from new and old networks
tfsNew = table(newNetTopEdges[,1])
tfsOld = table(araboxcisFiltered[,1])[names(tfsNew)]

# Histogram of degrees of TFs vs frequency
hist(as.numeric(tfsNew), 
     main = 'SinceAraBOXcis', 
     xlab = 'degree of TFs')

# Do the same TFs have high degrees in AraBOXcis and the new network? 
plot(as.numeric(tfsNew),
     as.numeric(tfsOld),
     xlab = 'Degree in SinceAraBOXcis',
     ylab = 'Degree in AraBOXcis')

# The 20 TFs with the highest degrees
sort(tfsNew, decreasing = TRUE) [1:20]


# Arabidopsis genes can be searched in arabidopsis.org - use the eFP browser tool. 

# Three different metrics of TF importance will be calculated below. 

install.packages('igraph')
library(igraph)

install.packages('network')
library(network)

# Load pheatmap, create a simple network
library(pheatmap)
simple_network <- graph_from_edgelist(as.matrix(newNetTopEdges[,c(1,2)]))



# Calculating node betweenness metric
node_betweenness_all <- betweenness(simple_network)
node_betweenness = node_betweenness_all[which(node_betweenness_all > 0)]
sort(node_betweenness, decreasing = TRUE) [1:20]

# Plots index vs node_betweenness sorted
plot(sort(node_betweenness))



# Calculating node centrality metric
node_centrality_all <- alpha_centrality(simple_network)
node_centrality = node_centrality_all[which(node_centrality_all > 0)]
sort(node_centrality, decreasing = TRUE) [1:20]

# Plots index vs node centrality metric
plot(sort(node_centrality))



# Calculating hub score metric
node_hub_all <- hub_score(simple_network)$vector
node_hub = node_hub_all[which(node_hub_all > 0)]
sort(node_hub, decreasing = TRUE) [1:20]

# Plots index vs node hub metric
sort(plot(node_hub))



# Comparing node betweenness vs node centrality
plot(node_betweenness_all, node_centrality_all)



# Comparing node hub vs node centrality
plot(node_hub_all, node_centrality_all)



