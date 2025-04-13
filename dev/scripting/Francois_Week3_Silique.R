#Silique-week 3
#1. LOADING NECESSARY DATA FROM WEEK1

#install.packages('Matrix')
#load libraries
library (tidyverse)
library(Matrix)
library(ggplot2)

#install packages
install.packages("BiocManager")

#Install devtools to load R helper functions
install.packages("devtools")
library(devtools)


file.exists('dev/utilities/dataprocessingHelperFunctions.R')
#install HelperFunctions.R
source('dev/utilities/dataprocessingHelperFunctions.R')

#load Silique data
a=load('data/GSE226097_silique_230221.RData')

#Load the original AraBOXcis network  trained on bulk RNA-seq in seedlings 
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

#Instaling and Loading packages for week3
BiocManager::install("GENIE3")


#SILIQUE-WEEK3

BiocManager::install("GENIE3")

library(GENIE3)
library(Matrix)

#load Silique data
a=load('data/GSE226097_silique_230221.RData')

#Load the original AraBOXcis network  trained on bulk RNA-seq in seedlings 
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

#View no of transcription factors
tfs=unique(araboxcis[,1])

#filter for transcription factors also present in scRNA seq data
tfSubs=tfs[which(tfs %in% rownames(gbox))]

net=GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees=5) #use more than 5 trees if your computer can handle it

save(net, file='Silique_network_nTree_5.RData')

net=load('Silique-network_nTree_5.RData')

load('Silique_network_nTree_5.RData')

str(net)

load('Silique_network_nTree_5.RData')


#Code from Github R Helper functions

convertToAdjacency <- function(network, threshold){
  
  temp=unlist(sapply(rownames(network), function(i){
    paste(i, colnames(network)[which(network[i,]>threshold)], network[i,which(network[i,]>threshold)])
  }))
  
  col3=data.frame(t(sapply(temp, function(i){
    split=strsplit(i, ' ')
    c(split[[1]][1], split[[1]][2], split[[1]][3])
  })))
  
  col3[,3]=as.numeric(col3[,3])
  
  col3
}

ginieOutput=convertToAdjacency(net, 0.05)


#Checking network dimensions
dim(ginieOutput)

#Checking first 10 edges in the network
ginieOutput[1:10,]

#loading the network
a=load('Silique_network_nTree_5.RData') ### Remember to actually load in your specific network 
newNet=GENIE3::getLinkList(net)
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

##Find overlaps in the Silique network and Arabidopsis network
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

##Find important genes in the network
tfsNew=table(newNetTopEdges[,1])
tfsOld=table(araboxcisFiltered[,1])[names(tfsNew)]

#histogram of degrees should look like an exponential distribution, because biological networks are a kind of network called 'scale-free', meaning that most TFs only regulate a small number of genes, but a few are influential hubs.
hist(as.numeric(tfsNew), main='SinceAraBOXcis', xlab='degree of TFs')

#Checking if the same TFs have high degrees in AraBOXcis and our new network:
plot(as.numeric(tfsNew), as.numeric(tfsOld), xlab='degree in SinceAraBOXcis', ylab='degree in AraBOXcis')

# print out the 20 TFs with highest degrees.  
sort(tfsNew, decreasing=TRUE)[1:20]

##calculate all the different metrics of TF importance

#loading more packages
library(igraph)
library(network)
library(pheatmap)

#Creating a Simple Network: constructing a graph  using graph_from_edgelist()
simple_network <- graph_from_edgelist(as.matrix(newNetTopEdges[,c(1,2)]))
node_betweenness_all <- betweenness(simple_network)
node_betweenness=node_betweenness_all[which(node_betweenness_all>0)]
sort(node_betweenness, decreasing=TRUE)[1:20]

#viewing centrality between Silique network and Arabidopsis network
plot(sort(node_betweenness))

# First, let's look at the size of your network
print(vcount(simple_network))  # Number of nodes
print(ecount(simple_network))  # Number of edges

# Instead of alpha_centrality, use other centrality measures
# Betweenness (which you were trying to analyze earlier)
betw <- betweenness(simple_network)

# Degree centrality (simpler to calculate)
deg <- degree(simple_network)

# Closeness centrality
clos <- closeness(simple_network)

# Eigenvector centrality (may be more stable than alpha)
try({
  eig <- eigen_centrality(simple_network)$vector
})

# Find top genes using degree (which should be very memory efficient)
sorted_deg <- sort(deg, decreasing=TRUE)
top_genes_deg <- names(sorted_deg[1:10])
print("Top 10 genes by degree:")
print(top_genes_deg)

# If betweenness worked, get top genes by betweenness
if(exists("betw")) {
  sorted_betw <- sort(betw, decreasing=TRUE)
  top_genes_betw <- names(sorted_betw[1:10])
  print("Top 10 genes by betweenness:")
  print(top_genes_betw)
}

# If you already have a list of TFs (tfSubs from your earlier code)
tf_list <- tfSubs

# Check which of your top genes are TFs
top_betweenness_tfs <- intersect(names(sorted_betw[1:10]), tf_list)
print("Top betweenness genes that are TFs:")
print(top_betweenness_tfs)





# Install Bioconductor if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install DOSE package from Bioconductor
BiocManager::install("DOSE")

#load DOSE
library(DOSE)

#Save betweenness,node-centrality,node-hub as centrality_Silique.Rdata
save(node_betweenness,node_centrality,node_hub,file='data/centrality_Silique.RData')

load(data/centrality_Silique.RData)
load("data/centrality_Silique.RData")
ls

node_betweenness
node_centrality
node_hub

summary(node_betweenness)
summary(node_centrality)
summary(node_hub)

summary(node_hub)

# Sort the betweenness centrality values in decreasing order and get the top 10
top_10_betweenness <- sort(node_betweenness, decreasing = TRUE)[1:10]

# Display the top 10 genes and their betweenness centrality values
top_10_betweenness

# Sort the node centrality values in decreasing order and get the top 10
top_10_centrality <- sort(node_centrality, decreasing = TRUE)[1:10]

# Display the top 10 genes and their node centrality values
top_10_centrality

# Sort the hub scores in decreasing order and get the top 10
top_10_hub <- sort(node_hub, decreasing = TRUE)[1:10]

# Display the top 10 genes and their hub scores
top_10_hub



#calculating and sorting hub scores in Silique network
#abline(h=5000)
node_hub_all <- hub_score(simple_network)$vector
node_hub=node_hub_all[which(node_hub_all>0)]
sort(node_hub, decreasing=TRUE)[1:20]

install.packages("igraph")  # Install igraph if necessary
library(igraph)  # Load the igraph package
library(igraph)
# Example edge list: Replace with your actual data
edges_df <- data.frame(
  from = c("AT3G07340", "AT4G29100", "AT5G24800", "AT3G59060", "AT3G54620"), # Gene A (source)
  to = c("AT5G24800", "AT3G59060", "AT3G54620", "AT3G24140", "AT5G49450")     # Gene B (target)
)

#view hub scores in Silique network
plot(sort(node_hub))

#abline(h=0.6)

plot(node_betweenness_all, node_centrality_all)

# Extract nodes that have values for both betweenness and alpha centrality
common_nodes <- intersect(names(node_betweenness_all), names(node_centrality_all))

# Subset both vectors to only include the common nodes
node_betweenness_common <- node_betweenness_all[common_nodes]
node_centrality_common <- node_centrality_all[common_nodes]

# Plot the values of betweenness vs. centrality
plot(node_betweenness_common, node_centrality_common,
     xlab = "node_hub_all", ylab = "node_centrality_all",
     main = "Betweenness vs.  Centrality")

# Get the names of the nodes from both vectors
common_nodes <- intersect(names(node_hub_all), names(node_centrality_all))

# Subset the vectors based on the common nodes
node_hub_common <- node_hub_all[common_nodes]
node_centrality_common <- node_centrality_all[common_nodes]

plot(node_hub_all, node_betweenness_all)

# 1. Extract top genes for each metric separately
# Assuming your network metrics are in dataframes or vectors named accordingly

# For centrality (top 20 genes)
top_centrality <- data.frame(gene = names(node_centrality_all), 
                             centrality = node_centrality_all)
top_centrality <- top_centrality[order(-top_centrality$centrality), ]
top_centrality_genes <- head(top_centrality$gene, 20)

# For centrality (top 20 genes)
top_centrality <- data.frame(gene = names(node_centrality_all), 
                             centrality = node_centrality_all)
top_centrality <- top_centrality[order(-top_centrality$centrality), ]
top_centrality_genes <- head(top_centrality$gene, 20)

# Then your intersection code should work
common_genes <- Reduce(intersect, list(top_centrality_genes, 
                                       top_betweenness_genes, 
                                       top_hub_genes))
library(biomaRt)
listMarts()

# Try with the correct URL protocol
ensembl_plants <- useMart(biomart="plants_mart",
                          host="http://plants.ensembl.org")

Good news! The connection seems to be working now, though there's a warning that Ensembl will soon require HTTPS. Let's proceed with your analysis:
  RCopy# Update to use HTTPS as recommended
ensembl_plants <- useMart(biomart="plants_mart",
                          host="https://plants.ensembl.org")

# Now let's list the available datasets
datasets <- listDatasets(ensembl_plants)

# Find the Arabidopsis dataset
arabidopsis_dataset <- datasets[grep("thaliana", datasets$dataset),]
print(arabidopsis_dataset)

# Connect to the Arabidopsis dataset
arabidopsis <- useDataset(arabidopsis_dataset$dataset[1], mart=ensembl_plants)

# Now you can query your genes
# First, let's check what attributes are available
attributes <- listAttributes(arabidopsis)
head(attributes)

# Check what filters are available
filters <- listFilters(arabidopsis)
head(filters)


BiocManager::biomaRt
library(biomaRt)
mart <- useMart("plants_mart", dataset="athaliana_eg_gene")
results <- getBM(attributes=c("ensembl_gene_id", "description", "go_id"), 
                 filters="ensembl_gene_id", 
                 values=your_gene_list, 
                 mart=mart)

##Finding associations between GO terms in Silique network
#load data on functional data required
a=load('data/functionalData.RData')

 #using HelperFunctions
 source('dev/utilities/dataprocessingHelperFunctions.R')
                                   
 #Using pafway from HelperFunction
 pafway <- function(GO, edges, GOtypes) {
 GOinNetwork = GO[unique(c(edges[, 1], edges[, 2]))]
grepLen=sapply(GOtypes, function(i){
 length(grep(i, GOinNetwork))
})
names(grepLen)=GOtypes
sapply(GOtypes, function(i) {
if(grepLen[i]!=0){

# find edges first:
grepI=grep(i, GOinNetwork[edges[, 1]])
sapply(GOtypes, function(j) {
if(grepLen[j]!=0){
a=length(grep(j, GOinNetwork[edges[grepI, 2]]))
p_bot = grepLen[i]/length(GOinNetwork) * grepLen[j]/length(GOinNetwork)
 b = stats::binom.test(a, length(edges[, 1]), p = p_bot, alternative = c("greater"))
b$p.value
}else{1}
})
}else{rep(1, length(GOtypes))}
})
                                   }
   
                                   
                                   #Running pathway function
                                   pafwayOut=pafway(GOconcat, newNetTopEdges, unique(goOfInterest))
                                   rownames(pafwayOut)=colnames(pafwayOut)
                                   
                                   #Filter to only include rows and columns with at least one significant factor:
                                   atLeastOneSigRow=which(apply(pafwayOut, 1, function(i){length(which(i<0.05))})>0)
                                   
                                   
                                   atLeastOneSigCol=which(apply(pafwayOut, 2, function(i){length(which(i<0.05))})>0)
                                   
                                   pafwayInterestingOnly=pafwayOut[atLeastOneSigRow, atLeastOneSigCol]
                                   
                                   #instaling packages for pheatmap
                                   install.packages("pheatmap")
                                   
                                   #loading pheatmap
                                   library(pheatmap)
                                   
                                   #Creating heatmap- columns are upstream of the terms associated with the rows.  
                                   #The values correspond to p-values, so smaller is more significant.
                                   pheatmap(pafwayInterestingOnly)
                                   
                                   #Creating another heatmap to zoom  most significant associations 
                                   #by taking a log
                                   pheatmap(log(pafwayInterestingOnly, 10))
                                   
                                   #improving labeling
                                   pheatmap(log(pafwayInterestingOnly, 10), 
                                            fontsize_row = 7)  
                                   
                                     
                                   