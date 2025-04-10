# CREATING GO HEATMAPS 
# AUTHOR: Zeynep Sakaoglu
# Part 1 --------------------------------------------------------------------------------------


# Loading necessary libraries
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(gridExtra)

# Defining all datasets
datasets <- c(
  "seedling_0d", 
  "seedling_3d", 
  "seedling_6d", 
  "seedling_12d", 
  "Rosette_21d", 
  "Rosette_30d", 
  "stem", 
  "Flower",
  "Silique"  
)

# Creating a function to load and process data
load_data <- function(dataset_name) {
  file_path <- paste0("exported_data/", dataset_name, "_go_results.csv")
  if (file.exists(file_path)) {
    data <- read.csv(file_path)
    data$Dataset <- dataset_name
    return(data)
  } else {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
}

# Load all datasets
all_go_data <- list()
for (dataset in datasets) {
  data <- load_data(dataset)
  if (!is.null(data)) {
    all_go_data[[dataset]] <- data
  }
}

# Checking for Na values and whether there is data in the data frame
if (length(all_go_data) == 0) {
  stop("No data was loaded. Please check file paths.")
}

# Combine all datasets
combined_go_results <- do.call(rbind, all_go_data)

# Get all unique GO terms
all_terms <- unique(c(
  unique(combined_go_results$DownstreamTerm),
  unique(combined_go_results$UpstreamTerm)
))

# Create a function to generate a matrix for each dataset
create_go_matrix <- function(go_data, all_terms) {

# Create an empty matrix
  matrix_data <- matrix(NA, nrow = length(all_terms), ncol = length(all_terms))
  rownames(matrix_data) <- all_terms
  colnames(matrix_data) <- all_terms
  
# Fill in the matrix with p-values
  for (i in 1:nrow(go_data)) {
    row_term <- go_data$DownstreamTerm[i]
    col_term <- go_data$UpstreamTerm[i]
    p_value <- go_data$pValue[i]
    
    row_idx <- which(all_terms == row_term)
    col_idx <- which(all_terms == col_term)
    
    if (length(row_idx) > 0 && length(col_idx) > 0) {
      matrix_data[row_idx, col_idx] <- p_value
    }
  }
  
  return(matrix_data)
}

# Create matrices for each dataset
dataset_matrices <- list()
for (dataset in names(all_go_data)) {
  dataset_matrices[[dataset]] <- create_go_matrix(
    filter(combined_go_results, Dataset == dataset),
    all_terms
  )
}


# Part 2 filtering for significance-------------------------------------------------------------------

# Filter for significant GO terms
significance_threshold <- 0.05

# Function to identify significant terms in a matrix
get_significant_terms <- function(matrix_data, threshold) {
  significant_rows <- apply(matrix_data, 1, function(x) any(x <= threshold, na.rm = TRUE))
  significant_cols <- apply(matrix_data, 2, function(x) any(x <= threshold, na.rm = TRUE))
  
  return(names(which(significant_rows | significant_cols)))
}

# Find significant terms in each dataset
sig_terms_list <- list()
for (dataset in names(dataset_matrices)) {
  sig_terms_list[[dataset]] <- get_significant_terms(dataset_matrices[[dataset]], significance_threshold)
  cat(paste0("Found ", length(sig_terms_list[[dataset]]), " significant GO terms in ", dataset, "\n"))
}

# Combine significant terms across all datasets
all_sig_terms <- unique(unlist(sig_terms_list))
cat(paste0("Total unique significant GO terms across all datasets: ", length(all_sig_terms), "\n"))

# 52 significant genes in total-


# Optimizing the content for visualization
max_terms <- 100  
if (length(all_sig_terms) > max_terms) {
# Calculate minimum p-value for each term across datasets
  min_pvalues <- numeric(length(all_sig_terms))
  names(min_pvalues) <- all_sig_terms
  
  for (i in 1:length(all_sig_terms)) {
    term <- all_sig_terms[i]
    p_values <- c()
    
    for (dataset in names(dataset_matrices)) {
      matrix_data <- dataset_matrices[[dataset]]
      
# Get p-values in both directions (as downstream and upstream)
      p_down <- min(matrix_data[term, ], na.rm = TRUE)
      p_up <- min(matrix_data[, term], na.rm = TRUE)
      
# Handle NA/Inf values
      p_down <- if(is.infinite(p_down)) 1 else p_down
      p_up <- if(is.infinite(p_up)) 1 else p_up
      
      p_values <- c(p_values, p_down, p_up)
    }
    
# Use minimum p-value across all datasets and positions
    min_pvalues[i] <- min(p_values)
  }
  
# Sort and select top terms
  all_sig_terms <- names(sort(min_pvalues)[1:max_terms])
  cat(paste0("Selected top ", max_terms, " significant GO terms for visualization\n"))
}

# Extract the sub-matrices with only significant terms
sig_matrices <- list()
for (dataset in names(dataset_matrices)) {
  
# Use intersect to ensure to only select terms that exist in the matrix
  available_terms <- intersect(all_sig_terms, rownames(dataset_matrices[[dataset]]))
  sig_matrices[[dataset]] <- dataset_matrices[[dataset]][available_terms, available_terms]
}

# Transform p-values to -log10(p) for better visualization (OPTIONAL/ BUT LOOKS BETTER THIS WAY)
transform_matrix <- function(matrix_data) {
  
# Replace NA with 1 (no significance) (OPTONAL)
  matrix_data[is.na(matrix_data)] <- 1
  
# Transform to -log10(p)
  log_matrix <- -log10(matrix_data)
  
# Cap values 
  max_value <- 10
  log_matrix[log_matrix > max_value] <- max_value
  
  return(log_matrix)
}

log_matrices <- list()
for (dataset in names(sig_matrices)) {
  log_matrices[[dataset]] <- transform_matrix(sig_matrices[[dataset]])
}


# Part 3-------------------------------------------------------------------------------------
# Create a combined matrix for all datasets
combined_log_matrix <- matrix(NA, 
                              nrow = length(all_sig_terms) * length(log_matrices), 
                              ncol = length(all_sig_terms))

all_dataset_names <- names(log_matrices)
row_names <- c()

# Fill the matrix and create row names
for (i in 1:length(all_dataset_names)) {
  dataset <- all_dataset_names[i]
  start_idx <- (i-1) * length(all_sig_terms) + 1
  end_idx <- i * length(all_sig_terms)
  
  combined_log_matrix[start_idx:end_idx, ] <- log_matrices[[dataset]]
  row_names <- c(row_names, paste0(dataset, "_", all_sig_terms))
}

rownames(combined_log_matrix) <- row_names
colnames(combined_log_matrix) <- all_sig_terms

# Create annotation for rows to indicate datasets
dataset_annotation <- data.frame(
  Dataset = rep(all_dataset_names, each = length(all_sig_terms))
)
rownames(dataset_annotation) <- row_names



# Part 4--------------------------------------------------------------------------------------
# Setup colors for annotation
num_datasets <- length(all_dataset_names)
dataset_colors <- brewer.pal(min(9, num_datasets), "Set1")
if (num_datasets > 9) {
  
# Adding colors part
  dataset_colors <- c(dataset_colors, brewer.pal(min(8, num_datasets-9), "Set2"))
}
names(dataset_colors) <- all_dataset_names[1:length(dataset_colors)]

annotation_colors <- list(
  Dataset = setNames(dataset_colors, all_dataset_names)
)

# Create individual heatmaps for each dataset for further use 
pdf("figures/individual_GO_heatmaps.pdf", width = 15, height = 15 * ceiling(length(log_matrices)/2))
par(mfrow = c(ceiling(length(log_matrices)/2), 2))

plot_list <- list()
for (dataset in names(log_matrices)) {
  p <- pheatmap(
    log_matrices[[dataset]],
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
    main = paste0("GO Enrichment: ", dataset),
    fontsize_row = 8,
    fontsize_col = 8,
    angle_col = 45,
    silent = TRUE
  )
  plot_list[[dataset]] <- p[[4]]
}
dev.off()

# Create the combined heatmap
pdf("figures/all_datasets_GO_heatmap.pdf", width = 20, height = 30)
pheatmap(
  combined_log_matrix,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6,  # Smaller font due to more rows
  fontsize_col = 8,
  angle_col = 45,
  annotation_row = dataset_annotation,
  annotation_colors = annotation_colors,
  main = "Gene Ontology Enrichment Across All Datasets",
  legend_labels = c("Not significant", "Highly significant"),
  legend_breaks = c(0, 10),
  legend = TRUE
)
dev.off()

# Arrange plots in a grid
n_plots <- length(plot_list)
n_cols <- 2
n_rows <- ceiling(n_plots/n_cols)

grid.arrange(grobs = plot_list, ncol = n_cols, nrow = n_rows)
dev.off()

# Creating aggregated maximum significance heatmap
# Create matrix for max value across all datasets
max_sig_matrix <- matrix(0, nrow = length(all_sig_terms), ncol = length(all_sig_terms))
rownames(max_sig_matrix) <- all_sig_terms
colnames(max_sig_matrix) <- all_sig_terms

for (dataset in names(log_matrices)) {
  max_sig_matrix <- pmax(max_sig_matrix, log_matrices[[dataset]])
}

pdf("figures/aggregated_max_GO_heatmap.pdf", width = 15, height = 15)
pheatmap(
  max_sig_matrix,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  main = "Maximum GO Enrichment Significance Across All Datasets",
  fontsize_row = 8,
  fontsize_col = 8,
  angle_col = 45,
  legend_breaks = c(0, 5, 10),
  legend_labels = c("Not significant", "Significant", "Highly significant")
)
dev.off()

# Create a combined similarity heatmap showing relationships between datasets
similarity_matrix <- matrix(0, 
                            nrow = length(log_matrices), 
                            ncol = length(log_matrices))
rownames(similarity_matrix) <- names(log_matrices)
colnames(similarity_matrix) <- names(log_matrices)

for (i in 1:length(names(log_matrices))) {
  for (j in 1:length(names(log_matrices))) {
    if (i != j) {
      
# Calculate correlation between datasets
      dataset1 <- log_matrices[[names(log_matrices)[i]]]
      dataset2 <- log_matrices[[names(log_matrices)[j]]]
      
# Flatten matrices for correlation
      flat1 <- as.vector(dataset1)
      flat2 <- as.vector(dataset2)
      
# Calculate correlation
      similarity_matrix[i, j] <- cor(flat1, flat2, method = "pearson")
    } else {
      similarity_matrix[i, j] <- 1  # Self-correlation
    }
  }
}

# FIX: Sort datasets by developmental stage order
# Create a vector for ordering
stage_values <- sapply(rownames(similarity_matrix), get_stage_order)
# Sort the rownames by their stage values
ordered_rownames <- rownames(similarity_matrix)[order(stage_values)]
# Now use these ordered names to subset the matrix
similarity_matrix <- similarity_matrix[ordered_rownames, ordered_rownames]

# Create the dataset similarity heatmap
pdf("figures/dataset_similarity_heatmap.pdf", width = 10, height = 8)
pheatmap(
  similarity_matrix,
  color = colorRampPalette(c("white", "darkblue"))(100),
  main = "Similarity Between Datasets (GO Term Enrichment Patterns)",
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 8
)
dev.off()

# Create a heatmap showing GO terms by dataset
go_by_dataset_matrix <- matrix(0, 
                               nrow = length(all_sig_terms), 
                               ncol = length(log_matrices))
rownames(go_by_dataset_matrix) <- all_sig_terms
colnames(go_by_dataset_matrix) <- names(log_matrices)

for (i in 1:length(all_sig_terms)) {
  term <- all_sig_terms[i]
  for (j in 1:length(names(log_matrices))) {
    dataset <- names(log_matrices)[j]
    
# Get maximum significance of this term in this dataset
    go_by_dataset_matrix[i, j] <- max(log_matrices[[dataset]][term, ], log_matrices[[dataset]][, term])
  }
}

# Order the columns by developmental stage
go_by_dataset_matrix <- go_by_dataset_matrix[, ordered_rownames]

# Create the GO terms by dataset heatmap
pdf("figures/go_terms_by_dataset_heatmap.pdf", width = 12, height = 20)
pheatmap(
  go_by_dataset_matrix,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  main = "GO Term Significance Across Developmental Stages",
  fontsize_row = 6,
  fontsize_col = 8,
  angle_col = 45,
  cluster_rows = TRUE,
  cluster_cols = FALSE  
)
dev.off()
# Keep developmental stage order- CLUSTER_COLS!


# Create heatmap with GO terms clustered by their patterns
pdf("figures/clustered_GO_terms_by_dataset.pdf", width = 12, height = 20)
pheatmap(
  go_by_dataset_matrix,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,  
  fontsize_row = 7,
  fontsize_col = 9,
  main = "GO Terms Clustered by Significance Pattern Across Development",
  annotation_col = data.frame(
    Stage = ordered_rownames,  
    row.names = ordered_rownames
  ),
  legend_breaks = c(0, 5, 10),
  legend_labels = c("Not significant", "Significant", "Highly significant")
)
dev.off()

# Create a network visualization of the most connected GO terms
# Find the most connected GO terms 
# Count connections with significance > 3

connectivity <- rowSums(max_sig_matrix > 3, na.rm = TRUE)  
top_connected <- names(sort(connectivity, decreasing = TRUE)[1:min(30, length(connectivity))])

# Extract sub-matrix for the most connected terms
connected_matrix <- max_sig_matrix[top_connected, top_connected]

# Save this as a heatmap
pdf("figures/highly_connected_GO_terms.pdf", width = 12, height = 10)
pheatmap(
  connected_matrix,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  main = "Most Connected GO Terms Across All Datasets",
  fontsize_row = 8,
  fontsize_col = 8,
  display_numbers = TRUE,
  number_format = "%.1f",
  fontsize_number = 7
)
dev.off()

# Print information about the results
cat("\n===== Gene Ontology Analysis Results =====\n")
cat(paste0("Analyzed ", length(log_matrices), " datasets\n"))
cat(paste0("Selected ", length(all_sig_terms), " significant GO terms\n"))
cat("\nFiles created in 'figures' directory:\n")
cat("1. all_datasets_GO_heatmap.pdf - Combined heatmap of all datasets\n")
cat("2. individual_GO_heatmaps.pdf - Individual heatmaps for each dataset\n")
cat("3. aggregated_max_GO_heatmap.pdf - Maximum significance across all datasets\n")
cat("4. GO_term_developmental_trajectory.pdf - Changes in top GO terms across developmental stages\n")
cat("5. dataset_similarity_heatmap.pdf - Similarity between datasets\n")
cat("6. clustered_GO_terms_by_dataset.pdf - GO terms clustered by significance patterns\n")
cat("7. highly_connected_GO_terms.pdf - Network of most connected GO terms\n")



# FIXING THE GRAPHS-----------------------------------------------------------------------------------------

# Fixing the GO_term_developmental_trajectory Plot
# Fixing the GO Term trajectory plot showing NA on the x-axis

# Recreate trajectory data from scratch with explicit dataset names
trajectory_data <- data.frame()

# Selecting top 20 most significant GO terms
top_terms <- names(sort(apply(max_sig_matrix, 1, max), decreasing = TRUE)[1:20])

# Explicitly define the ordered datasets
ordered_datasets <- c(
  "seedling_0d", 
  "seedling_3d", 
  "seedling_6d", 
  "seedling_12d", 
  "Rosette_21d", 
  "Rosette_30d", 
  "stem", 
  "Flower",
  "Silique"  
)

# Verify these datasets exist in log_matrices
available_datasets <- names(log_matrices)
print("Available datasets in log_matrices:")
print(available_datasets)

# Find any missing datasets
missing_datasets <- setdiff(ordered_datasets, available_datasets)
if(length(missing_datasets) > 0) {
  print("Warning: Some expected datasets are missing:")
  print(missing_datasets)
}

# Use only available datasets
ordered_datasets <- intersect(ordered_datasets, available_datasets)

# Create the trajectory data with explicit x-axis ordering
for (term in top_terms) {
  for (dataset in ordered_datasets) {
    
# Skip if dataset doesn't exist in log_matrices
    if (!dataset %in% names(log_matrices)) {
      next
    }
    
# Check if term exists in this dataset's matrix
    if (!term %in% rownames(log_matrices[[dataset]]) || 
        !term %in% colnames(log_matrices[[dataset]])) {
      
# If term doesn't exist, use 0 (not significant)
      sig_value <- 0
    } else {
      
# Get maximum significance for this term (either as upstream or downstream)
      sig_value <- max(log_matrices[[dataset]][term, ], 
                       log_matrices[[dataset]][, term], 
                       na.rm = TRUE)
      
# Handle Inf values
      if (is.infinite(sig_value)) sig_value <- 0
    }
    
# Add to trajectory data
    trajectory_data <- rbind(trajectory_data, data.frame(
      Term = term,
      Dataset = dataset,
      Significance = sig_value,
      stringsAsFactors = FALSE
    ))
  }
}

# Set the Dataset column as a factor with explicit ordering
trajectory_data$Dataset <- factor(trajectory_data$Dataset, 
                                  levels = ordered_datasets)

# Check the structure of our new data frame
print("Structure of new trajectory_data:")
str(trajectory_data)
print("Summary of Dataset factor levels:")
summary(trajectory_data$Dataset)

# Create the fixed trajectory plot with guaranteed x-axis labels
pdf("figures/fixed_GO_trajectory.pdf", width = 14, height = 10)

# Create the plot with explicit axis settings
p <- ggplot(trajectory_data, aes(x = Dataset, y = Significance, group = Term, color = Term)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    title = "GO Term Significance Across Developmental Stages",
    x = "Developmental Stage",
    y = "-log10(p-value)",
    color = "GO Term"
  ) +
  scale_y_continuous(limits = c(0, 10)) +
  scale_x_discrete(drop = FALSE)

print(p)
dev.off()

#-----------------------------------------------------------------------------------------------

# Adding x-axis jitter to the GO term trajectory plot

# Using the trajectory_data we created in the previous fix
# If you need to recreate it, include the data creation code here

# Create a trajectory plot with jitter
pdf("figures/GO_trajectory_with_jitter.pdf", width = 14, height = 10)

# Option 1: Using position_jitter for the points only
p1 <- ggplot(trajectory_data, aes(x = Dataset, y = Significance, group = Term, color = Term)) +
  # Keep lines without jitter to show the trend
  geom_line(size = 0.7, alpha = 0.7) +
  # Add jittered points
  geom_point(
    size = 3, 
    position = position_jitter(width = 0.2, height = 0, seed = 123)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    title = "GO Term Significance Across Developmental Stages (with Jitter)",
    x = "Developmental Stage",
    y = "-log10(p-value)",
    color = "GO Term"
  ) +
  scale_y_continuous(limits = c(0, 10)) +
  scale_x_discrete(drop = FALSE)

print(p1)
dev.off()


# Fixing graphs:
# Modified pheatmap call with additional parameters
pdf("figures/all_datasets_GO_heatmap.pdf", width = 20, height = 30)
pheatmap(
  combined_log_matrix,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6,  
  fontsize_col = 8,
  angle_col = 45,
  annotation_row = dataset_annotation,
  annotation_colors = annotation_colors,
  main = "Gene Ontology Enrichment Across All Datasets",
  legend_labels = c("Not significant", "Highly significant"),
  legend_breaks = c(0, 10),
  legend = TRUE,
  na_col = "grey90",  
  border_color = NA,  
  cutree_rows = 5,    
  cutree_cols = 5,
  treeheight_row = 100,
  treeheight_col = 100,
  cellwidth = NA,
  cellheight = NA,
  scale = "none",
  drop_levels = TRUE
)
dev.off()


# Fixing the trajectory and enrichment plots:
# Explicitly define ordered datasets
ordered_datasets <- c(
  "seedling_0d", 
  "seedling_3d", 
  "seedling_6d", 
  "seedling_12d", 
  "Rosette_21d", 
  "Rosette_30d", 
  "stem", 
  "Flower",
  "Silique"  
)

# Create trajectory data with explicit x-axis ordering
trajectory_data <- data.frame()

for (term in top_terms) {
  for (dataset in ordered_datasets) {
    # Safe handling of dataset and term existence
    if (!dataset %in% names(log_matrices)) next
    
    if (!term %in% rownames(log_matrices[[dataset]]) || 
        !term %in% colnames(log_matrices[[dataset]])) {
      sig_value <- 0
    } else {
      sig_value <- max(log_matrices[[dataset]][term, ], 
                       log_matrices[[dataset]][, term], 
                       na.rm = TRUE)
    }
    
    trajectory_data <- rbind(trajectory_data, data.frame(
      Term = term,
      Dataset = dataset,
      Significance = sig_value,
      stringsAsFactors = FALSE
    ))
  }
}

# Set Dataset as a factor with explicit ordering
trajectory_data$Dataset <- factor(trajectory_data$Dataset, 
                                  levels = ordered_datasets)

p1 <- ggplot(trajectory_data, aes(x = Dataset, y = Significance, group = Term, color = Term)) +
  geom_line(size = 0.7, alpha = 0.7) +
  geom_point(
    size = 3, 
    position = position_jitter(width = 0.2, height = 0, seed = 123)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    title = "GO Term Significance Across Developmental Stages (with Jitter)",
    x = "Developmental Stage",
    y = "-log10(p-value)",
    color = "GO Term"
  ) +
  scale_y_continuous(limits = c(0, 10)) +
  scale_x_discrete(drop = FALSE)  # Ensures all x-axis levels are shown



# final

# Create a trajectory plot with jitter and larger text
pdf("figures/GO_trajectory_with_jitter.pdf", width = 16, height = 12)

p1 <- ggplot(trajectory_data, aes(x = Dataset, y = Significance, group = Term, color = Term)) +
  # Keep lines without jitter to show the trend
  geom_line(size = 0.7, alpha = 0.7) +
  # Add jittered points
  geom_point(
    size = 3, 
    position = position_jitter(width = 0.2, height = 0, seed = 123)
  ) +
  theme_bw() +
  theme(
    # Increase text sizes
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Larger x-axis labels
    axis.text.y = element_text(size = 12),  # Larger y-axis labels
    axis.title.x = element_text(size = 14, face = "bold"),  # Larger x-axis title
    axis.title.y = element_text(size = 14, face = "bold"),  # Larger y-axis title
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Larger plot title
    legend.title = element_text(size = 12, face = "bold"),  # Larger legend title
    legend.text = element_text(size = 10),  # Larger legend text
    legend.position = "right"
  ) +
  labs(
    title = "GO Term Significance Across Developmental Stages",
    x = "Developmental Stage",
    y = "-log10(p-value)",
    color = "GO Term"
  ) +
  scale_y_continuous(limits = c(0, 10)) +
  scale_x_discrete(drop = FALSE)

print(p1)
dev.off()


# Verify these datasets exist in log_matrices
available_datasets <- names(log_matrices)
print("Available datasets in log_matrices:")
print(available_datasets)

# Find any missing datasets
missing_datasets <- setdiff(ordered_datasets, available_datasets)
if(length(missing_datasets) > 0) {
  print("Warning: Some expected datasets are missing:")
  print(missing_datasets)
}

# Use only available datasets
ordered_datasets <- intersect(ordered_datasets, available_datasets)

# Create the trajectory data with explicit x-axis ordering
for (term in top_terms) {
  for (dataset in ordered_datasets) {
    
    # Skip if dataset doesn't exist in log_matrices
    if (!dataset %in% names(log_matrices)) {
      next
    }
    
    # Check if term exists in this dataset's matrix
    if (!term %in% rownames(log_matrices[[dataset]]) || 
        !term %in% colnames(log_matrices[[dataset]])) {
      
      # If term doesn't exist, use 0 (not significant)
      sig_value <- 0
    } else {
      
      # Get maximum significance for this term (either as upstream or downstream)
      sig_value <- max(log_matrices[[dataset]][term, ], 
                       log_matrices[[dataset]][, term], 
                       na.rm = TRUE)
      
      # Handle Inf values
      if (is.infinite(sig_value)) sig_value <- 0
    }
    
    # Add to trajectory data
    trajectory_data <- rbind(trajectory_data, data.frame(
      Term = term,
      Dataset = dataset,
      Significance = sig_value,
      stringsAsFactors = FALSE
    ))
  }
}

# Set the Dataset column as a factor with explicit ordering
trajectory_data$Dataset <- factor(trajectory_data$Dataset, 
                                  levels = ordered_datasets)

# Check the structure of our new data frame
print("Structure of new trajectory_data:")
str(trajectory_data)
print("Summary of Dataset factor levels:")
summary(trajectory_data$Dataset)

# Create the fixed trajectory plot with guaranteed x-axis labels
pdf("figures/fixed_GO_trajectory.pdf", width = 14, height = 10)

# Create the plot with explicit axis settings
p <- ggplot(trajectory_data, aes(x = Dataset, y = Significance, group = Term, color = Term)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    title = "GO Term Significance Across Developmental Stages",
    x = "Developmental Stage",
    y = "-log10(p-value)",
    color = "GO Term"
  ) +
  scale_y_continuous(limits = c(0, 10)) +
  scale_x_discrete(drop = FALSE)

print(p)
dev.off()

#-----------------------------------------------------------------------------------------------

# Adding x-axis jitter to the GO term trajectory plot

# Using the trajectory_data we created in the previous fix
# If you need to recreate it, include the data creation code here

# Create a trajectory plot with jitter
pdf("figures/GO_trajectory_with_jitter.pdf", width = 14, height = 10)

# Option 1: Using position_jitter for the points only
p1 <- ggplot(trajectory_data, aes(x = Dataset, y = Significance, group = Term, color = Term)) +
  # Keep lines without jitter to show the trend
  geom_line(size = 0.7, alpha = 0.7) +
  # Add jittered points
  geom_point(
    size = 3, 
    position = position_jitter(width = 0.2, height = 0, seed = 123)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    title = "GO Term Significance Across Developmental Stages (with Jitter)",
    x = "Developmental Stage",
    y = "-log10(p-value)",
    color = "GO Term"
  ) +
  scale_y_continuous(limits = c(0, 10)) +
  scale_x_discrete(drop = FALSE)

print(p1)
dev.off()


# Modified pheatmap of all_datasets call with additional parameters
pdf("figures/all_datasets_GO_heatmap.pdf", width = 20, height = 30)
pheatmap(
  combined_log_matrix,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6,  
  fontsize_col = 8,
  angle_col = 45,
  annotation_row = dataset_annotation,
  annotation_colors = annotation_colors,
  main = "Gene Ontology Enrichment Across All Datasets",
  legend_labels = c("Not significant", "Highly significant"),
  legend_breaks = c(0, 10),
  legend = TRUE,
  na_col = "grey90",  
  border_color = NA,  
  cutree_rows = 5,    
  cutree_cols = 5,
  treeheight_row = 100,
  treeheight_col = 100,
  cellwidth = NA,
  cellheight = NA,
  scale = "none",
  drop_levels = TRUE
)
dev.off()