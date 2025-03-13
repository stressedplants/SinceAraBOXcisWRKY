# Week 3 Homework:  

# MAKING A NETWORK FROM SCRATCH

BiocManager::install("GENIE3")
library(GENIE3)

# Helper Function
source("C:/Users/User/OneDrive/Documents/Seed_0days_Cookeygam/utilities.R")

# Loading single cell data sets ('Seed_0days')
Seed0Days = load("C:/Users/User/OneDrive/Documents/Seed_0days_Cookeygam/GSE226097_seed_0d_230221.RData")

# Loading the original AraBoxcis network that was trained on bulk RNA-seq in seedlings
AraBoxcis = read.csv("C:/Users/User/OneDrive/Documents/Seed_0days_Cookeygam/gboxNetwork22C.csv", header=TRUE)

# A list of Transcription factors
Tfs = unique(AraBoxcis[,1])

# A list of the regulators (TfSubs)
TfSubs = Tfs[which(Tfs %in% rownames(gbox))]

net = GENIE3(as.matrix(gbox),
             regulators = TfSubs, nTrees = 5)
save(net, file = 'Seed0Days')

str(net)

ginieOutput = convertToAdjacency(net, 0.05)

# Dimensions of first 10 edges network
dim(ginieOutput)
ginieOutput[1: 10,]

# Loading the Network

NewNet = GENIE3::getLinkList(net)

# Find Overlaps between the SingleCellGene & AraBoxcis network
# Gets the set of unique genes in your new netwrok
GenesInNet = unique(c(NewNet[,1], NewNet[,2]))

# Filter AraBoxcis Network only to genes in new network
 AraBoxcisFiltered = AraBoxcis[which(AraBoxcis[,1] %in% GenesInNet & AraBoxcis[,2] %in% GenesInNet),]
 
 # Extract the top edges in new network to make network same size
 NewNetTopEdges = NewNet[1: length(AraBoxcisFiltered[,1]),]
 
 # Reformat edges so it is more straightforward to compare them
 EdgesNew = paste(NewNetTopEdges[,1], NewNetTopEdges[,2], sep = '_')
 
 head(EdgesNew)
 EdgesOld = paste(AraBoxcisFiltered[,1], AraBoxcisFiltered[,2], sep = '_')
 # Overlap
 length(which(EdgesNew %in% EdgesOld))

 # New network only
 length(which(! (EdgesNew %in% EdgesOld)))
 
 # Old Newtwork only
 length(! (which(EdgesOld %in% EdgesNew)))
 
 
 # Making a ven diagram
 # Installing matplot
 install.packages('VennDiagram')
 library(VennDiagram)
 
New_network_only = 41048
Old_network_only = 8497
Overlap_edges =  8497

# creating a Venn Diagram
venn.plot = draw.pairwise.venn(
  area1 = New_network_only + Overlap_edges,
  area2 = Old_network_only + Overlap_edges,
  cross.area = Overlap_edges,
  category = c('New Network(Seed 0 Days)', 'Old Network(AraBoxcis)'),
  fill = c("#F28E8E", "#98C798"),
  alpha = 0.7,
  lty = 'solid',
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(-20, 20)
)


grid.draw(venn.plot)

# Important genes
TfsNew = table(NewNetTopEdges[,1])
TfsOld = table(AraBoxcisFiltered[,1])[names(TfsNew)]

# VIRTUALIZATIONS
# Histogram
hist(as.numeric(TfsNew), col = 'skyblue', main = 'SinceAraboxcis', xlab = 'degree of TFs')

# scatter plot
# Define the custom colors
custom_colors = c("#FF00FF", "turquoise4", "seagreen", "goldenrod3")

# Convert tfsNew and tfsOld to numeric
TfsNew_numeric = as.numeric(TfsNew)
TfsOld_numeric = as.numeric(TfsOld)

# groups based on quartiles of tfsNew
degree_groups = cut(TfsNew_numeric, breaks = 4, labels = custom_colors)

# Plot with color-coded groups
plot(TfsNew_numeric, TfsOld_numeric, 
     xlab = 'degree in SinceAraBOXcis', 
     ylab = 'degree in AraBOXcis', 
     col = as.character(degree_groups), 
   
     main = "Gene Network Degree Comparison")

# Adding a legend
legend("topright", legend = c("Low", "Medium-Low", "Medium-High", "High"), 
       col = custom_colors, pch = 16, title = "Degree Groups")

# Print out the 20 TFs with highest degrees. 
sort(TfsNew, decreasing = TRUE)[1:20]

# installing packages for calculating 
# all the different metrics of TF importance
install.packages('igraph')
library(igraph)


install.packages('network')
library(network)


install.packages('pheatmap')
library('pheatmap')

simple_network = graph_from_edgelist(as.matrix(NewNetTopEdges[,c(1,2)]))

node_betweenness_all = betweenness(simple_network)
node_betweenness = node_betweenness_all[which(node_betweenness_all > 0)]
sort(node_betweenness_all, decreasing = TRUE)[1:20]

# Plotting Nodes "BETWEENESS"
plot(sort(node_betweenness), col = 'green',)



# abline(h = 5000)

node_centrality_all = alpha_centrality(simple_network, alpha=0.9)
node_centrality=node_centrality_all[which(node_centrality_all>0)]
sort(node_centrality, decreasing=TRUE)[1:20]


# NODE CENTRALITY
plot(sort(node_centrality), col = 'blue')

node_hub_all = hub_score(simple_network)$vector
node_hub = node_hub_all[which(node_hub_all>0)]     
sort(node_hub, decreasing = TRUE)[1:20]

plot(sort(node_hub), col = 'purple')

plot(node_betweenness_all, node_centrality_all)
length(node_betweenness_all)
length(node_centrality_all)
str(node_centrality_all)
summary(node_centrality_all)


plot(node_hub_all, node_centrality_all, col = "red" )

# Association between GO terms in the network

a = load("C:/Users/User/OneDrive/Documents/Seed_0days_Cookeygam/functionalData - Copy.RData")
source("C:/Users/User/OneDrive/Documents/Seed_0days_Cookeygam/utilities.R")

pafwayOut = pafway(GOconcat, NewNetTopEdges, unique(goOfInterest))
rownames(pafwayOut) = colnames(pafwayOut)

# Filter to only include row and columns with at least one significant factor:
atleastOneSigRow = which(apply(pafwayOut, 1, function(i){length(which(i<0.05))})>0)

atleastOneSigCol = which(apply(pafwayOut, 2, function(i){length(which(i<0.05))})>0)

pafwayIntrestingOnly = pafwayOut[atleastOneSigRow, atleastOneSigCol]

#Here, the terms associated with the columns are upstream of the terms associated with the rows. 
#The values correspond to p-values, so smaller is more significant.

install.packages("pheatmap")
library('pheatmap')
pheatmap(pafwayIntrestingOnly)

#Let's re-do zooming to the most significant associations by taking a log

pheatmap(log(pafwayIntrestingOnly, 10))
