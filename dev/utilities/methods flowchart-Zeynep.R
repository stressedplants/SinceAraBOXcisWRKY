# Install if not already installed
install.packages("DiagrammeR")

# Load package
library(DiagrammeR)

# Create the flowchart
DiagrammeR::grViz("
digraph G {
  
  graph [layout = dot, rankdir = TB]

  # Define nodes
  node [shape = rectangle, style = filled, fontname = Arial]
  
  DataAcq [label = 'Data Acquisition\nLoad scRNA-seq & AraBOXcis data', fillcolor = LightBlue]
  QC [label = 'Quality Control\nFilter low-quality cells and genes', fillcolor = LightBlue]
  DimRed [label = 'Dimensionality Reduction\nUMAP & PCA visualization', fillcolor = LightBlue]
  TFSelect [label = 'TF Selection\nFilter TFs present in both datasets', fillcolor = LightBlue]
  NetInf [label = 'Network Inference\nGENIE3 (nTrees = 5, 8, 10)', fillcolor = LightCoral]
  NetComp [label = 'Network Comparison\nVenn diagrams & edge overlap', fillcolor = LightGreen]
  CentAnal [label = 'Centrality Analysis\nBetweenness, Alpha, Hub scores', fillcolor = LightGreen]
  PathAnal [label = 'Pathway Analysis\nGO terms & p-value heatmaps', fillcolor = LightGreen]
  IntInterp [label = 'Integration & Interpretation\nKey regulator identification', fillcolor = Orange]
  
  # Define edges
  DataAcq -> QC
  QC -> DimRed
  DimRed -> TFSelect
  TFSelect -> NetInf
  NetInf -> NetComp
  NetInf -> CentAnal
  NetInf -> PathAnal
  NetComp -> IntInterp
  CentAnal -> IntInterp
  PathAnal -> IntInterp
}
")

# Install if needed
install.packages(c("ggplot2", "ggplotify"))

# Load libraries
library(ggplot2)
library(ggplotify)

#=============================


# Install necessary packages if not installed
install.packages(c("ggplot2", "ggrepel", "gridExtra"))

# Load libraries
library(ggplot2)
library(ggrepel)

# Define flowchart elements with better positioning
flowchart_data <- data.frame(
  x = c(3, 3, 3, 3, 3, 1.5, 3, 4.5, 3),
  y = c(9, 8, 7, 6, 5, 4, 4, 4, 3),
  text = c("Data Acquisition", "Quality Control", "Dimensionality Reduction", 
           "TF Selection", "Network Inference", "Network Comparison", 
           "Centrality Analysis", "Pathway Analysis", "Integration & Interpretation"),
  color = c("blue", "blue", "blue", "blue", "red", "green", "green", "green", "orange")
)

# Create flowchart
ggplot(flowchart_data, aes(x, y)) +
  geom_tile(aes(fill = color), width = 1.5, height = 0.5, color = "black") +
  geom_text(aes(label = text), size = 5, fontface = "bold", color = "black") +
  scale_fill_manual(values = c("blue" = "lightblue", "red" = "lightcoral",
                               "green" = "lightgreen", "orange" = "orange")) +
  theme_void() +
  theme(legend.position = "none") +
  # Adding arrows for better structure
  geom_segment(aes(x = 3, xend = 3, y = 8.7, yend = 8.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 3, xend = 3, y = 7.7, yend = 7.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 3, xend = 3, y = 6.7, yend = 6.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 3, xend = 3, y = 5.7, yend = 5.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 3, xend = 1.9, y = 4.7, yend = 4.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 3, xend = 3, y = 4.7, yend = 4.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 3, xend = 4.1, y = 4.7, yend = 4.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 1.5, xend = 3, y = 3.7, yend = 3.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 3, xend = 3, y = 3.7, yend = 3.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 4.5, xend = 3, y = 3.7, yend = 3.3), arrow = arrow(type = "closed"), size = 1)



# with color coding labels 
#-----------------------------

# Define flowchart elements with improved positioning
flowchart_data <- data.frame(
  x = c(3, 3, 3, 3, 3, 1.5, 3, 4.5, 3),
  y = c(9, 8, 7, 6, 5, 4, 4, 4, 3),
  text = c("Data Acquisition", "Quality Control", "Dimensionality Reduction", 
           "TF Selection", "Network Inference", "Network Comparison", 
           "Centrality Analysis", "Pathway Analysis", "Integration & Interpretation"),
  color = c("blue", "blue", "blue", "blue", "red", "green", "green", "green", "orange"),
  stringsAsFactors = FALSE
)

# Create flowchart
ggplot(flowchart_data, aes(x, y)) +
  geom_tile(aes(fill = color), width = 1.5, height = 0.5, color = "black") +
  geom_text(aes(label = text), size = 4, fontface = "bold", color = "black") +
  scale_fill_manual(values = c("blue" = "lightblue", "red" = "lightcoral",
                               "green" = "lightgreen", "orange" = "orange")) +
  theme_void() +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  scale_color_manual(values = c("blue" = "lightblue", "red" = "lightcoral", 
                                "green" = "lightgreen", "orange" = "orange")) +
  # Adding arrows for better structure
  geom_segment(aes(x = 3, xend = 3, y = 8.7, yend = 8.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 3, xend = 3, y = 7.7, yend = 7.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 3, xend = 3, y = 6.7, yend = 6.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 3, xend = 3, y = 5.7, yend = 5.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 3, xend = 1.9, y = 4.7, yend = 4.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 3, xend = 3, y = 4.7, yend = 4.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 3, xend = 4.1, y = 4.7, yend = 4.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 1.5, xend = 3, y = 3.7, yend = 3.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 3, xend = 3, y = 3.7, yend = 3.3), arrow = arrow(type = "closed"), size = 1) + 
  geom_segment(aes(x = 4.5, xend = 3, y = 3.7, yend = 3.3), arrow = arrow(type = "closed"), size = 1) +
  # Add legend for the color coding
  scale_fill_manual(values = c("blue" = "lightblue", "red" = "lightcoral", "green" = "lightgreen", "orange" = "orange"),
                    labels = c("Data Processing", "Network Generation", "Network Analysis", "Integration"))




