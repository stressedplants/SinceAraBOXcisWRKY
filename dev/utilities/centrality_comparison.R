# =============================================================================
# Title:        Centrality Comparison Across Developmental Stages
# Description:  This R script compares different types of network centrality 
#               measures—Alpha Centrality (AC), Betweenness Centrality (BC), and 
#               Degree Centrality (DC)—across various developmental stages in 
#               Arabidopsis. The script loads precomputed centrality metrics, 
#               compiles them into data frames, and generates log-log plots to 
#               visualize correlations between centrality measures. Custom 
#               themes and color palettes are included for standard and 
#               poster-ready outputs.
#
# Author:       Gagan  
# Date:         24-03-2025
#
# Citation:     If you use this script, please cite:
#               Gagan, "Centrality Comparison Across Developmental Stages", 2025.
# =============================================================================

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

