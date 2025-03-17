# Data import and organisation --------------------------------------------

# create an object with file path
data_list <- c('data/centrality/centrality_seedlingd12.RData',
               'data/centrality/centrality_seedlingd6.RData',
               'data/centrality/centrality_flower.RData',
               'data/centrality/centrality_rosette21D.RData',
               'data/centrality/centrality_Silique.RData',
               'data/centrality/centrality_rosette30.RData',
               'data/centrality/centrality_seedd0.RData')

# create an object with names
name_list <- c('sdd12','sdd6','flowr','rsd21','silqu','rsd30','seed')

# initialise empty lists
acs <- list()
bcs <- list()
dcs <- list()



for (i in seq_along(data_list)) { # loop over the lenght of data_list
  # load the data
  load(data_list[i])
  
  # generate the variable names
  ac_name <- paste0("ac_", name_list[i])
  bc_name <- paste0("bc_", name_list[i])
  dc_name <- paste0("dc_", name_list[i])
  
  # assign the loaded data to generated variable name
  assign(ac_name, node_centrality)
  assign(bc_name, node_betweenness)
  assign(dc_name, node_hub)
  
  #create a list of variable names
  acs <- append(acs, ac_name)
  bcs <- append(bcs, bc_name)
  dcs <- append(dcs, dc_name)
  
}

remove(ac_name,bc_name,dc_name,node_betweenness,node_centrality,node_hub, 
       i,data_list,name_list) #delete unnecessary variables


# Creating a combined table -------------------------------------------------

library(tidyverse)

# top n genes to be looked at
n=10


# for alpha centrality
acs_top <- c()

for (x in acs) {
  x <- sort(get(x),decreasing = T)[1:n]
  acs_top <- unique(c(acs_top,names(x)))
}

acs_top <- data_frame(Gene = acs_top)

for (x in acs) {
  y <- get(x)
  acs_top <- acs_top %>% 
    left_join(data.frame(Gene = names(y), x=y), by = "Gene") %>% 
    rename_with(~ x, .cols = x)
}

acs_top <- as.data.frame(acs_top)
row.names(acs_top) <- acs_top$Gene
acs_top$Gene <- NULL


# for betweenness centrality
bcs_top <- c()

for (x in bcs) {
  x <- sort(get(x),decreasing = T)[1:n]
  bcs_top <- unique(c(bcs_top,names(x)))
}

bcs_top <- data_frame(Gene = bcs_top)

for (x in bcs) {
  y <- get(x)
  bcs_top <- bcs_top %>% 
    left_join(data.frame(Gene = names(y), x=y), by = "Gene") %>% 
    rename_with(~ x, .cols = x)
}

bcs_top <- as.data.frame(bcs_top)
row.names(bcs_top) <- bcs_top$Gene
bcs_top$Gene <- NULL


# for degree centrality

dcs_top <- c()

for (x in dcs) {
  x <- sort(get(x),decreasing = T)[1:n]
  dcs_top <- unique(c(dcs_top,names(x)))
}

dcs_top <- data_frame(Gene = dcs_top)

for (x in dcs) {
  y <- get(x)
  dcs_top <- dcs_top %>% 
    left_join(data.frame(Gene = names(y), x=y), by = "Gene") %>% 
    rename_with(~ x, .cols = x)
}

dcs_top <- as.data.frame(dcs_top)
row.names(dcs_top) <- dcs_top$Gene
dcs_top$Gene <- NULL

remove(n,x,y)

# Presence Table  ---------------------------------------------------------

# create a truth table
acs_top_pres <- lapply(acs_top, function(x) as.numeric(!is.na(x))) %>% as.data.frame()
rownames(acs_top_pres) <- row.names(acs_top)
bcs_top_pres <- lapply(bcs_top, function(x) as.numeric(!is.na(x))) %>% as.data.frame()
rownames(bcs_top_pres) <- row.names(bcs_top)
dcs_top_pres <- lapply(dcs_top, function(x) as.numeric(!is.na(x))) %>% as.data.frame()
rownames(dcs_top_pres) <- row.names(dcs_top)

# plot a heatmap from the truth table
library(pheatmap)

pheatmap(acs_top_pres,main = 'Alpha Centrality'
         ,cluster_rows = FALSE, cluster_cols = FALSE, 
         color = c("white", "green"))

pheatmap(bcs_top_pres,main = 'Betweenness Centrality',
         cluster_rows = FALSE, cluster_cols = FALSE, 
         color = c("white", "green"))

pheatmap(dcs_top_pres,main = 'Degree Centrality',
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = c("white", "green"))

# genes present in all
acs_pres_all <- rownames(acs_top_pres)[rowSums(acs_top_pres == 1) == ncol(acs_top_pres)]
bcs_pres_all <- rownames(bcs_top_pres)[rowSums(bcs_top_pres == 1) == ncol(bcs_top_pres)]
dcs_pres_all <- rownames(dcs_top_pres)[rowSums(dcs_top_pres == 1) == ncol(dcs_top_pres)]


# Heat map ----------------------------------------------------------------

# change NA values to 0
acs_top_new <- acs_top
acs_top_new[is.na(acs_top_new)] <- 0

bcs_top_new <- bcs_top
bcs_top_new[is.na(bcs_top_new)] <- 0

dcs_top_new <- dcs_top
dcs_top_new[is.na(dcs_top_new)] <- 0

#plot heat maps
pheatmap(acs_top_new, main = 'Alpha Centrality'
         ,cluster_rows = FALSE, cluster_cols = FALSE)

pheatmap(bcs_top_new, main = 'Betweenness Centrality'
         ,cluster_rows = FALSE, cluster_cols = FALSE)

pheatmap(dcs_top_new, main = 'Degree Centrality'
         ,cluster_rows = FALSE, cluster_cols = FALSE)


# create log values
acs_top_new_log <- log1p(acs_top_new)

bcs_top_new_log <- log1p(bcs_top_new)

dcs_top_new_log <- log1p(dcs_top_new)

#plot log heat maps
pheatmap(acs_top_new_log, main = 'Log Alpha Centrality'
         ,cluster_rows = FALSE, cluster_cols = FALSE)

pheatmap(bcs_top_new_log, main = 'Log Betweenness Centrality'
         ,cluster_rows = FALSE, cluster_cols = FALSE)

pheatmap(dcs_top_new_log, main = 'Log Degree Centrality'
         ,cluster_rows = FALSE, cluster_cols = FALSE)


# Output Files ------------------------------------------------------------

write.csv(acs_top_pres, file = 'data/centrality/output/acs_pres.csv')
write.csv(bcs_top_pres, file = 'data/centrality/output/bcs_pres.csv')
write.csv(dcs_top_pres, file = 'data/centrality/output/dcs_pres.csv')

write_lines(acs_pres_all,'data/centrality/output/acs_pres_all.txt')
write_lines(bcs_pres_all,'data/centrality/output/bcs_pres_all.txt')
write_lines(dcs_pres_all,'data/centrality/output/dcs_pres_all.txt')

# write(paste(acs_pres_all, collapse = ","),
#       "data/centrality/output/acs_pres_all.txt")
# write(paste(bcs_pres_all, collapse = ","),
#       "data/centrality/output/bcs_pres_all.txt")
# write(paste(dcs_pres_all, collapse = ","),
#       "data/centrality/output/dcs_pres_all.txt")