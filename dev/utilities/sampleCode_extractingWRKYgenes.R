library(Seurat)
#read in full file 
a=LoadSeuratRds("data/GSE226097_flower_230221.rds")

#read in WRKY TFs
#source: Eulgem, T.,Rushton, P. J.,Robatzek, S.,Somssich, I. E. (2000) The WRKY superfamily of plant transcription factors. TRENDS IN PLANT SCIENCE 5
tfs=read.csv("data/WRKY_TFs.csv", header=F)[,1]

#filter for only TFs with over 100 transcripts found
length(which(rowSums(a[["RNA"]])[tfs]>100))
tfsInData=tfs[which(rowSums(a[["RNA"]])[tfs]>100)]

#read in genes near W-boxes
#source: https://doi.org/10.1038/s41598-019-38757-7
genesNearWbox=read.csv("data/WBoxgenes.csv", header=T)

#filter for only genes with over 100 transcripts found
length(which(rowSums(a[["RNA"]])[genesNearWbox[,1]]>100))
genesInData=genesNearWbox[,1][which(rowSums(a[["RNA"]])[genesNearWbox[,1]]>100)]

#draw histograms of read counts
hist(rowSums(a[["RNA"]])[tfsInData])
hist(rowSums(a[["RNA"]])[genesInData])

#extract only WRKY related data
incl=colSums(a[["RNA"]]$data[tfsInData,])
subSet=a[["RNA"]]$data[unique(c(tfsInData, genesInData)), which(incl>0)]

write.table(as.matrix(subSet), file="data/flowerSubset_Wbox.csv")
