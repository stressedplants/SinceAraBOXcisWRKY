# =============================================================================
# Title:        GO Term Enrichment Visualization for Network Centrality-Based Gene Sets
# Description:  This R script analyzes and visualizes the distribution of GO terms
#               associated with genes ranked by Alpha Centrality, Betweenness Centrality,
#               and Degree Centrality. It parses TAIR-derived annotations and a custom
#               functional dataset to identify and plot the most frequent biological 
#               processes within each high-centrality gene group.
#
#               The output includes bar plots of top GO terms for each centrality class.
#
# Author:       Gagan  
# Date:         24-03-2025
#
# Citation:     If you use this script, please cite:
#               Gagan, "GO Term Enrichment Visualization for Network 
#               Centrality-Based Gene Sets", 2025.
# =============================================================================

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


