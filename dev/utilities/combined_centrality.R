# Data import and organisation --------------------------------------------

# create an object with file path
data_list <- c('data/centrality/centrality_seedd0.RData',
               'data/centrality/centrality_seedlingd3.RData',
               'data/centrality/centrality_seedlingd6.RData',
               'data/centrality/centrality_seedlingd12.RData',
               'data/centrality/centrality_rosette21D.RData',
               'data/centrality/centrality_rosette30.RData',
               'data/centrality/centrality_stem.RData',
               'data/centrality/centrality_Silique.RData',
               'data/centrality/centrality_flower.RData')

# create an object with names
name_list <- c('seed','sdd3','sdd6','sdd12','rsd21','rsd30','stem','silqu','flowr')

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

pheatmap(acs_top_pres,main = 'Alpha Centrality',
         color = c("white", "brown"))

pheatmap(bcs_top_pres,main = 'Betweenness Centrality',
         color = c("white", "brown"))

pheatmap(dcs_top_pres,main = 'Degree Centrality',
         color = c("white", "brown"))

# genes present in all
acs_pres_all <- rownames(acs_top_pres)[rowSums(acs_top_pres == 1) == ncol(acs_top_pres)]
bcs_pres_all <- rownames(bcs_top_pres)[rowSums(bcs_top_pres == 1) == ncol(bcs_top_pres)]
dcs_pres_all <- rownames(dcs_top_pres)[rowSums(dcs_top_pres == 1) == ncol(dcs_top_pres)]

# genes present in more than 6 stages
acs_rowsum <- rowSums(acs_top_pres)
bcs_rowsum <- rowSums(bcs_top_pres)
dcs_rowsum <- rowSums(dcs_top_pres)

acs_pres_maj <- rownames(acs_top_pres)[which(acs_rowsum>6)]
bcs_pres_maj <- rownames(bcs_top_pres)[which(bcs_rowsum>6)]
dcs_pres_maj <- rownames(dcs_top_pres)[which(dcs_rowsum>6)]


# Heat map ---rowSums()# Heat map ----------------------------------------------------------------

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

write_lines(acs_pres_maj, 'data/centrality/output/acs_pres_maj.txt')
write_lines(bcs_pres_maj, 'data/centrality/output/bcs_pres_maj.txt')
write_lines(dcs_pres_maj, 'data/centrality/output/dcs_pres_maj.txt')


# Centrality comparison ---------------------------------------------------

## Loading the new data ----------------------------------------------------

# create an object with file path
data_list <- c('data/centrality/All/centrality_all_seedlingd12.RData',
               'data/centrality/All/centrality_all_seedling0days.RData',
               'data/centrality/All/centrality_all_Rosette21D.RData',
               'data/centrality/All/centrality_all_seedlingd6.RData',
               'data/centrality/All/centrality_all_rosette30.RData',
               'data/centrality/All/centrality_all_seedlingd3.RData',
               'data/centrality/All/centrality_all_Silique.RData')

# create an object with names
name_list <- c('sdd12','seed','rsd21','sdd6','rsd30','sdd3','silqu')

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
  assign(ac_name, node_centrality_all)
  assign(bc_name, node_betweenness_all)
  assign(dc_name, node_hub_all)
  
  #create a list of variable names
  acs <- append(acs, ac_name)
  bcs <- append(bcs, bc_name)
  dcs <- append(dcs, dc_name)
  
}

remove(ac_name,bc_name,dc_name,node_betweenness_all,node_centrality_all,node_hub_all, 
       i,data_list,name_list) #delete unnecessary variables

## AC vs BC ----------------------------------------------------------------

# Load necessary package
library(ggplot2)

# Create a data frame with all points
ACvBC <- data.frame(
  x = c(ac_sdd12, ac_sdd6, ac_seed, ac_rsd21, ac_rsd30,ac_sdd3,ac_silqu),
  y = c(bc_sdd12, bc_sdd6, bc_seed, bc_rsd21, bc_rsd30,bc_sdd3,bc_silqu),
  group = rep(c("Seedling Day12", "Seedling Day6", "Seed", "Rosette Day21", "Rosette Day30","Seedling Day3", "Silique"),
              times = c(length(ac_sdd12), length(ac_sdd6), length(ac_seed), length(ac_rsd21), length(ac_rsd30),length(ac_sdd3),length(ac_silqu)))
)

# Plot using ggplot2 with log-log scaling
ggplot(ACvBC, aes(x = x, y = y, color = group, shape = group)) +
  geom_point() +
  scale_x_log10() +  # Log scale for x-axis
  scale_y_log10() +  # Log scale for y-axis
  theme_minimal() +
  labs(title = 'Alpha Centrality Vs Betweenness Centrality',
       x = "Alpha Centrality Score (log scale)", y = "Betweenness Centrality Score (log scale)",
       color = "Stage", shape = "Stage") +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange","chocolate4", "plum3"))+  # Custom colors
  scale_shape_manual(values = c(16, 17, 18, 15, 8, 9, 11))  # Custom shapes

## DC v BC -----------------------------------------------------------------
###
# Create a data frame with all points
DCvBC <- data.frame(
  x = c(dc_sdd12, dc_sdd6, dc_seed, dc_rsd21, dc_rsd30, dc_sdd3,dc_silqu),
  y = c(bc_sdd12, bc_sdd6, bc_seed, bc_rsd21, bc_rsd30, bc_sdd3,bc_silqu),
  group = rep(c("Seedling Day12", "Seedling Day6", "Seed", "Rosette Day21", "Rosette Day30", "Seedling Day3", "Silique"),
              times = c(length(dc_sdd12), length(dc_sdd6), length(dc_seed), length(dc_rsd21), length(dc_rsd30), length(dc_sdd3),length(dc_silqu))
))

# Plot using ggplot2 with log-log scaling
ggplot(DCvBC, aes(x = x, y = y, color = group,shape = group)) +
  geom_point(size = 1.5) +
  scale_x_log10(limits = c(0.00001, 1)) +  # Log scale for x-axis
  scale_y_log10() +  # Log scale for y-axis
  theme_minimal() +
  labs(title = 'Degree Centrality Vs Betweenness Centrality',
       x = "Degree Centrality Score (log scale)", y = "Betweenness Centrality Score (log scale)",
       color = "Stage", shape = "Stage") +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange","chocolate4","plum3"))+  # Custom colors
  scale_shape_manual(values = c(16, 17, 18, 15, 8, 9, 11))  # Custom shapes
""