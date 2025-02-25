

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



