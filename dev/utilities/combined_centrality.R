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


## poster recolour ---------------------------------------------------------
pheatmap(acs_top_pres,main = 'Alpha Centrality',
         color = c("white", "brown"))
pheatmap(bcs_top_pres,main = 'Betweenness Centrality',
         color = c("white", "brown"))



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
               'data/centrality/All/centrality_all_Silique.RData',
               'data/centrality/All/centrality_all_flower.RData',
               'data/centrality/All/centrality_all_stem.RData')

# create an object with names
name_list <- c('sdd12','seed','rsd21','sdd6','rsd30','sdd3','silqu','flowr','stem')

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
  x = c(ac_sdd12, ac_sdd6, ac_seed, ac_rsd21, ac_rsd30,ac_sdd3,ac_silqu,ac_flowr,ac_stem),
  y = c(bc_sdd12, bc_sdd6, bc_seed, bc_rsd21, bc_rsd30,bc_sdd3,bc_silqu,bc_flowr,bc_stem),
  group = rep(c("Seedling Day12", "Seedling Day6", "Seed", "Rosette Day21", 
                "Rosette Day30","Seedling Day3", "Silique","Flower","Stem"),
              times = c(length(ac_sdd12), length(ac_sdd6), length(ac_seed),
                        length(ac_rsd21), length(ac_rsd30),length(ac_sdd3),
                         length(ac_silqu),length(ac_flowr),length(ac_stem)))
)

# Plot using ggplot2 with log-log scaling
ggplot(ACvBC, aes(x = x, y = y, color = group, shape = group)) +
  geom_point(size = 2.5) +
  scale_x_log10() +  # Log scale for x-axis
  scale_y_log10() +  # Log scale for y-axis
  theme_minimal() +
  labs(title = 'Alpha Centrality Vs Betweenness Centrality',
       x = "Alpha Centrality Score (log scale)", 
       y = "Betweenness Centrality Score (log scale)",
       color = "Stage", shape = "Stage") +
  scale_color_manual(values = c("red", "blue", "green", "purple",
                                "orange","chocolate4", "plum3",
                                "orangered","cyan3"))+  # Custom colors
  scale_shape_manual(values = c(16, 17, 18, 15, 8, 9, 11,7,13))+  # Custom shapes
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=25),
        plot.title=element_text(size=40),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20))


### Poster recolour ---------------------------------------------------------
# Plot recolour
ggplot(ACvBC, aes(x = x, y = y, color = group, shape = group)) +
  geom_point(size = 2.5) +
  scale_x_log10() +  # Log scale for x-axis
  scale_y_log10() +  # Log scale for y-axis
  theme_minimal() +
  labs(title = 'Alpha Centrality Vs Betweenness Centrality',
       x = "Alpha Centrality Score (log scale)", 
       y = "Betweenness Centrality Score (log scale)",
       color = "Stage", shape = "Stage") +
  scale_color_manual(values = c("red", "blue", "green", "purple",
                                "orange","chocolate4", "plum3",
                                "orangered","cyan3"))+  # Custom colors
  scale_shape_manual(values = c(16, 17, 18, 15, 8, 9, 11,7,13))+  # Custom shapes
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=25),
        plot.title=element_text(size=40),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20))+
  theme(plot.background = element_rect(fill = '#759D5D'),
        panel.grid.major=element_line(colour="#3a4e2e"),
        panel.grid.minor=element_line(colour="#3a4e2e"),
        plot.title = element_text(colour = "#ffffff"),
        axis.title = element_text(colour = "#ffffff"),
        axis.text = element_text(colour = "#ffffff"),
        legend.text=element_text(colour = "#ffffff"),
        legend.title=element_text(colour = "#ffffff"))


## DC v BC -----------------------------------------------------------------
# Create a data frame with all points
DCvBC <- data.frame(
  x = c(dc_sdd12, dc_sdd6, dc_seed, dc_rsd21, dc_rsd30, dc_sdd3,dc_silqu,dc_flowr,dc_stem),
  y = c(bc_sdd12, bc_sdd6, bc_seed, bc_rsd21, bc_rsd30, bc_sdd3,bc_silqu,bc_flowr,bc_stem),
  group = rep(c("Seedling Day12", "Seedling Day6", "Seed", "Rosette Day21", 
                "Rosette Day30", "Seedling Day3", "Silique","Flower","Stem"),
              times = c(length(dc_sdd12), length(dc_sdd6), length(dc_seed), 
                        length(dc_rsd21), length(dc_rsd30), length(dc_sdd3),
                        length(dc_silqu),length(dc_flowr),length(dc_stem))
))

# Plot using ggplot2 with log-log scaling
ggplot(DCvBC, aes(x = x, y = y, color = group,shape = group)) +
  geom_point(size = 2.5) +
  scale_x_log10(limits = c(0.00001, 1)) +  # Log scale for x-axis
  scale_y_log10() +  # Log scale for y-axis
  theme_minimal() +
  labs(title = 'Degree Centrality Vs Betweenness Centrality',
       x = "Degree Centrality Score (log scale)", 
       y = "Betweenness Centrality Score (log scale)",
       color = "Stage", shape = "Stage") +
  scale_color_manual(values = c("red", "blue", "green", "purple",
                                "orange", "chocolate4", "plum3",
                                "orangered", "cyan3"))+  # Custom colors
  scale_shape_manual(values = c(16, 17, 18, 15, 8, 9, 11,7,13))+  # Custom shapes
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=25),
        plot.title=element_text(size=40),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20))


### Poster recolour ---------------------------------------------------------

# Plot recolour
ggplot(DCvBC, aes(x = x, y = y, color = group,shape = group)) +
  geom_point(size = 2.5) +
  scale_x_log10(limits = c(0.00001, 1)) +  # Log scale for x-axis
  scale_y_log10() +  # Log scale for y-axis
  theme_minimal() +
  labs(title = 'Degree Centrality Vs Betweenness Centrality',
       x = "Degree Centrality Score (log scale)", 
       y = "Betweenness Centrality Score (log scale)",
       color = "Stage", shape = "Stage") +
  scale_color_manual(values = c("red", "blue", "green", "purple",
                                "orange", "chocolate4", "plum3",
                                "orangered", "cyan3"))+  # Custom colors
  scale_shape_manual(values = c(16, 17, 18, 15, 8, 9, 11,7,13))+  # Custom shapes
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=25),
        plot.title=element_text(size=40),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20))+
  theme(plot.background = element_rect(fill = '#759D5D'),
        panel.grid.major=element_line(colour="#3a4e2e"),
        panel.grid.minor=element_line(colour="#3a4e2e"),
        plot.title = element_text(colour = "#ffffff"),
        axis.title = element_text(colour = "#ffffff"),
        axis.text = element_text(colour = "#ffffff"),
        legend.text=element_text(colour = "#ffffff"),
        legend.title=element_text(colour = "#ffffff"))


# GO terms analysis using TAIR--------------------------------------------------
#source('dev/utilities/dataprocessingHelperFunctions.R')

library(readr)
library(tidyverse)

num = 25

## ACs ---------------------------------------------------------------------
acs_maj <- read_tsv('data/centrality/TAIR_search_res/ACs_maj.tsv')
acs_maj <- acs_maj[,c("Locus","Keywords")]

acs_maj_kc<- acs_maj %>%
  separate_rows(Keywords, sep = ";") %>%  # Split on ';' delimiter
  count(Keywords, sort = TRUE)  # Count occurrences

# Select top num keywords
acs_maj_kc_top <- acs_maj_kc %>% 
  slice_max(n, n = num)  # Get top 10 most frequent keywords

# Plot frequency graph
ggplot(acs_maj_kc_top, aes(x = reorder(Keywords, n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip axes for better readability
  labs(title = "High Alpha Centrality Genes", x = "GO Terms", y = "Frequency") +
  theme_minimal()



## BCs ---------------------------------------------------------------------
bcs_maj <- read_tsv('data/centrality/TAIR_search_res/BCs_maj.tsv')
bcs_maj <- bcs_maj[,c("Locus","Keywords")]

bcs_maj_kc<- bcs_maj %>%
  separate_rows(Keywords, sep = ";") %>%  # Split on ';' delimiter
  count(Keywords, sort = TRUE)  # Count occurrences

# Select top num keywords
bcs_maj_kc_top <- bcs_maj_kc %>% 
  slice_max(n, n = num)  # Get top 10 most frequent keywords

# Plot frequency graph
ggplot(bcs_maj_kc_top, aes(x = reorder(Keywords, n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip axes for better readability
  labs(title = "High Betweenness Genes", x = "GO Terms", y = "Frequency") +
  theme_minimal()


## DCs ---------------------------------------------------------------------

dcs_maj <- read_tsv('data/centrality/TAIR_search_res/DCs_maj.tsv')
dcs_maj <- dcs_maj[,c("Locus","Keywords")]

dcs_maj_kc<- dcs_maj %>%
  separate_rows(Keywords, sep = ";") %>%  # Split on ';' delimiter
  count(Keywords, sort = TRUE)  # Count occurrences

# Select top num keywords
dcs_maj_kc_top <- dcs_maj_kc %>% 
  slice_max(n, n = num)  # Get top 10 most frequent keywords

# Plot frequency graph
ggplot(dcs_maj_kc_top, aes(x = reorder(Keywords, n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip axes for better readability
  labs(title = "High Degree Genes", x = "GO Terms", y = "Frequency") +
  theme_minimal()


# GO term analysis using given functional data ----------------------------

library(readr)
library(tidyverse)
library(stringr)


load('data/functionalData.RData')

GO_df <- data.frame(
  Locus = names(GOconcat), 
  Keywords = unlist(GOconcat), 
  stringsAsFactors = FALSE
)

num = 5

## ACs ---------------------------------------------------------------------

acs_maj <- read_lines('data/centrality/output/acs_pres_maj.txt')

# Filter for Locus values present in acs_maj
acs_maj_fn <- GO_df %>% 
  filter(Locus %in% acs_maj)

# Count occurrences of goOfInterest terms in acs_maj_fn$Keywords
acs_gooi_cnt <- sapply(goOfInterest, function(term) {
  sum(str_detect(acs_maj_fn$Keywords, fixed(term, ignore_case = TRUE)))
})


# Convert to a data frame
acs_gooi_cnt <- data.frame(Term = names(acs_gooi_cnt), Frequency = acs_gooi_cnt)
acs_gooi_cnt <- acs_gooi_cnt %>% arrange(desc(Frequency))

# Plot top num terms
ggplot(acs_gooi_cnt %>% slice_max(Frequency, n = num), 
       aes(x = reorder(Term, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top GO terms in genes with high Alpha Centrality",
       x = "GO Terms",
       y = "Frequency of GO term") +
  theme_minimal()


## BCs ---------------------------------------------------------------------

bcs_maj <- read_lines('data/centrality/output/bcs_pres_maj.txt')

# Filter for Locus values present in bcs_maj
bcs_maj_fn <- GO_df %>% 
  filter(Locus %in% bcs_maj)

# Count occurrences of goOfInterest terms in bcs_maj_fn$Keywords
bcs_gooi_cnt <- sapply(goOfInterest, function(term) {
  sum(str_detect(bcs_maj_fn$Keywords, fixed(term, ignore_case = TRUE)))
})


# Convert to a data frame
bcs_gooi_cnt <- data.frame(Term = names(bcs_gooi_cnt), Frequency = bcs_gooi_cnt)
bcs_gooi_cnt <- bcs_gooi_cnt %>% arrange(desc(Frequency))

# Plot top num terms
ggplot(bcs_gooi_cnt %>% slice_max(Frequency, n = num), 
       aes(x = reorder(Term, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "brown") +
  coord_flip() +
  labs(title = "Top GO terms in genes with high Betweenness Centrality",
       x = "GO Terms",
       y = "Frequency of GO term") +
  theme_minimal()+
  scale_y_continuous(limits = c(0,10))+
  scale_x_discrete(labels = function(x) tools::toTitleCase(x))+
  theme(axis.text=element_text(size=20,
                               margin = margin(, b = 30, l = 50)),
        axis.title=element_text(size=25),
        plot.title=element_text(size=35),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20))


### Poster recolour ---------------------------------------------------------

# Replot with colours
ggplot(bcs_gooi_cnt %>% slice_max(Frequency, n = num), 
       aes(x = reorder(Term, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "#275B42") +
  coord_flip() +
  labs(title = "Top GO terms in genes with high Betweenness Centrality",
       x = "GO Terms",
       y = "Frequency of GO term") +
  theme_minimal()+
  scale_y_continuous(limits = c(0,10))+
  scale_x_discrete(labels = function(x) tools::toTitleCase(x))+
  theme(axis.text=element_text(size=20,
                               margin = margin(, b = 30, l = 50)),
        axis.title=element_text(size=25),
        plot.title=element_text(size=35),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20))+
  theme(plot.background = element_rect(fill = '#759D5D'),
        panel.grid.major=element_line(colour="#3a4e2e"),
        panel.grid.minor=element_line(colour="#3a4e2e"),
        plot.title = element_text(colour = "#ffffff"),
        axis.title = element_text(colour = "#ffffff"),
        axis.text = element_text(colour = "#ffffff"))


## DCs ---------------------------------------------------------------------

dcs_maj <- read_lines('data/centrality/output/dcs_pres_maj.txt')

# Filter for Locus values present in dcs_maj
dcs_maj_fn <- GO_df %>% 
  filter(Locus %in% dcs_maj)

# Count occurrences of goOfInterest terms in dcs_maj_fn$Keywords
dcs_gooi_cnt <- sapply(goOfInterest, function(term) {
  sum(str_detect(dcs_maj_fn$Keywords, fixed(term, ignore_case = TRUE)))
})


# Convert to a data frame
dcs_gooi_cnt <- data.frame(Term = names(dcs_gooi_cnt), Frequency = dcs_gooi_cnt)
dcs_gooi_cnt <- dcs_gooi_cnt %>% arrange(desc(Frequency))

# Plot top num terms
ggplot(dcs_gooi_cnt %>% slice_max(Frequency, n = num), 
       aes(x = reorder(Term, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "brown") +
  coord_flip() +
  labs(title = "Top GO terms in genes with high Degree Centrality",
       x = "GO Terms",
       y = "Frequency of GO term") +
  theme_minimal()+
  scale_x_discrete(labels = function(x) tools::toTitleCase(x))+
  theme(axis.text=element_text(size=20,
                               margin = margin(, b = 30, l = 50)),
        axis.title=element_text(size=25),
        plot.title=element_text(size=35),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20))

### Poster recolour ---------------------------------------------------------

# Replot with colours
ggplot(dcs_gooi_cnt %>% slice_max(Frequency, n = num), 
       aes(x = reorder(Term, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "#275B42") +
  coord_flip() +
  labs(title = "Top GO terms in genes with high Degree Centrality",
       x = "GO Terms",
       y = "Frequency of GO term") +
  theme_minimal()+
  scale_x_discrete(labels = function(x) tools::toTitleCase(x))+
  theme(axis.text=element_text(size=20,
                               margin = margin(, b = 30, l = 50)),
        axis.title=element_text(size=25),
        plot.title=element_text(size=35),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),)+
  theme(plot.background = element_rect(fill = '#759D5D'),
        panel.grid.major=element_line(colour="#3a4e2e"),
        panel.grid.minor=element_line(colour="#3a4e2e"),
        plot.title = element_text(colour = "#ffffff"),
        axis.title = element_text(colour = "#ffffff"),
        axis.text = element_text(colour = "#ffffff"))


