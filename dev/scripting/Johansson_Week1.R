
# install.packages('Matrix')
library(Matrix)

# Loading helper functions
source('dev/utilities/dataprocessingHelperFunctions.R')

# if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

# Loading the flower data
a = load("data/flowerdata.RData", verbose=TRUE)

# Loading the original AraBOXcis network trained on bulk RNA-seq in seedlings
araboxcis <- read.csv(file = 'data/gboxNetwork22C.csv',header = TRUE)



# Contents of flower data
print(a)

# Dimensions of gbox
dim(gbox)

# Length of clust
length(clust)

# Names of the first 10 rows of gbox
rownames(gbox)[1:10]

# gbox column names
colnames(gbox)[1:10]

# First 50 rows and first 3 cols as matrix
as.matrix(gbox[1:50, 1:3])



# Cluster designations of the first 8 cells
clust[1:8]

# Cell counts in each cluster
table(clust)

# Plotting cell count against the clusters
plot(table(sort(clust)), 
     xlab = 'Cluster Name', 
     ylab = 'Number of Cells', 
     main = 'Histogram for the Flower Data')



# Dimensions of the Araboxcis matrix
dim(araboxcis)

# First 4 rows of the Araboxcis matrix
araboxcis[1:4,]

# Plotting a histogram of the edge scores
hist(araboxcis[,3])



# A list of the Transcription Factors 
tfs = unique(araboxcis[,1])

# Filtered list for only TFs also included in the single-cell RNA-seq dataset
tfSubs = tfs[which(tfs %in% rownames(gbox))]

# Number of TFs
length(tfSubs)



# Filter Genes and Cells With Very Low Values 

# Dimensions of the gbox
dim(gbox)

# Filtering-out cells with <1% of the genes expressed at all
thresh = 0.01
numberGenesPerCell = apply(gbox, 
                           2, 
                           function(i) {
                             length(which(i>0))
                             }
                           )
includeCells = which(numberGenesPerCell > (thresh * dim(gbox)[1])) 

gbox_filtered = gbox[,includeCells] 

dim(gbox_filtered)

# Removing genes that are expressed in <1% of the cells
numberCellPerGene = apply(gbox, 
                          1, 
                          function(i) {
                            length(which(i>0))
                            }
                          )
includeGenes = which(numberCellPerGene > (thresh * dim(gbox)[2]))

gbox_filtered = gbox_filtered[includeGenes,]

dim(gbox_filtered)



# First Visualisation With UMAP

install.packages('umap')
library(umap)

# Visualising data using UMAP
gbox.umap <- umap(gbox_filtered)

# Cell clustering of G-box related genes
colours = rainbow(length(unique(clust)))

plot(gbox.umap$layout[,1], 
     gbox.umap$layout[,2],
     col = colours[clust[includeCells]], 
     pch = 20,
     main = 'UMAP Flower', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2')


# PCA prior to UMAP:

pca = prcomp(gbox_filtered, 
             scale. = T, 
             rank. = 5)
gbox.pca.umap <- umap(pca$x)

colours = rainbow(length(unique(clust)))

plot(gbox.pca.umap$layout[,1], 
     gbox.pca.umap$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20, 
     main = 'PCA UMAP Flower', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2')



