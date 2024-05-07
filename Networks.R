# This is a file for diversity plots as well as network analysis

# Loading the packages
library(NetCoMi)
library(phyloseq)
library(igraph)
library(limma)
library(LaplacesDemon)
library(dplyr)
library(microbiome)
library(ggplot2)
library(ggforce)
library(sjPlot)
library(vegan)

# Path to the working directory
working_path = "/Users/apple/Desktop/Diversity/Network/"

# Loading the files and creation of phyloseq objects
environments <- c("all","Hospital", "MetaSUB", "Office")

for (env in environments) {
  # Construct file paths
  otu_file <- paste(working_path, env, 
                    "_filtered_data.csv", sep="")
  metadata_file <- paste(working_path, env, 
                         "_filtered_metadata.csv", sep="")
  taxa_file <- paste(working_path, env, 
                     "_filtered_taxa.csv", sep="")
  
  # Read OTU data
  otu_data <- read.table(otu_file, sep = ',', header=TRUE, 
                         check.names=FALSE)
  otu_data <- data.frame(otu_data[,-1], row.names = otu_data[,1], 
                         check.names=FALSE)
  otu_data <- otu_table(otu_data, taxa_are_rows = TRUE)
  
  # Read metadata
  meta_data <- read.table(metadata_file, sep = ',', header=TRUE)
  meta_data <- data.frame(meta_data[,-1], row.names = meta_data[,1])
  meta_data <- sample_data(meta_data)
  
  # Read taxonomy data
  taxa_data <- read.table(taxa_file, sep = ',', header=TRUE)
  taxa_data <- data.frame(taxa_data[,-1], row.names = taxa_data[,1])
  taxa_data <- as.matrix(taxa_data)
  taxa_data <- tax_table(taxa_data)
  
  physeq <- phyloseq(otu_data, taxa_data, meta_data) 
  assign(paste0("physeq_", env), physeq)
}

# Diversity analysis

# Alpha diversity
alpha_div <- microbiome::alpha(physeq_all, index = "all")

# Plot alpha diversity
p.shannon <- boxplot_alpha(physeq_all, 
                           index = "shannon",
                           x_var = "Environment")

p.shannon <- p.shannon + theme_minimal() + 
  labs(x="Environment", y="Shannon diversity") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))
p.shannon

# Testing differences
dif <- meta(physeq_all)
dif$diversity <- alpha_div$diversity_shannon
spl <- split(dif$diversity, dif$Environment)

p_matrix <- matrix(NA, nrow = length(spl), ncol = length(spl))

for (i in 1:length(spl)) {
  for (j in 1:length(spl)) {
    if (i != j) {
      ks_result <- ks.test(spl[[i]], spl[[j]])$p.value
      ks_result <- p.adjust(ks_result, method = "BH")
      p_matrix[i, j] <- ks_result
    } else {
      p_matrix[i, j] <- 1  # Diagonal elements (self-comparisons)
    }
  }
}

rownames(p_matrix) <- colnames(p_matrix) <- names(spl)

# Saving the plot
ggsave(paste0(working_path, "shannon_plot.svg"), plot = p.shannon, 
       width = 4, height = 5)

# Beta diversity
beta_div <- distance(physeq_all, method = "bray")

set.seed(123)
pcoa_result <- ordinate(physeq_all, method = "NMDS", distance = "bray")

p.beta <- plot_ordination(physeq_all, pcoa_result, color = "Environment") + 
  geom_point(aes(color = Environment), size = 3) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, group = Environment),
               geom = "polygon", alpha = 0.1, 
               linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 16)
  )

p.beta

ggsave(paste0(working_path, "beta_plot.svg"), plot = p.beta, 
              width = 6, height = 5)

# Testing differences
adonis_result <- adonis2(beta_div ~ Environment, 
                         data = as(sample_data(physeq_all), "data.frame"),
                         permutations = 999, method = "bray")


# Class-level plot
pseq.class.rel <-aggregate_rare(physeq_all, level = "class", 
                                detection = 0.001, prevalence = 0.1)


p.compositions <- plot_composition(pseq.class.rel,
                      average_by = "Environment") + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = "Study", 
       y = "Relative abundance")
#title = "Relative abundance data", 
#subtitle = "Subtitle",
#caption = "Caption text.") 
p.compositions <- p.compositions + scale_fill_brewer("Class", 
                      palette = "Paired") + theme_minimal()

p.compositions

ggsave(paste0(working_path, "composition_plot.svg"), plot = p.compositions, 
       width = 6, height = 5)

# Network construction
environments <- environments[environments != "all"]

for (env in environments) {
  print(env)
  physeq_env <- get(paste0("physeq_", env))
  
  # Network construction
  micro_net <- netConstruct(physeq_env,
                            measure = "gcoda",
                            normMethod = "clr",
                            zeroMethod = "multRepl",
                            dissFunc = "signed",
                            sparsMethod = "bootstrap",
                            alpha = 0.05,
                            adjust = "adaptBH",
                            nboot = 1000L,
                            assoBoot = NULL,
                            logFile = "log.txt",
                            cores = 3,
                            knnMutual = FALSE,
                            seed = 123,
                            verbose = 3)
  
  net_prop <- netAnalyze(micro_net, clustMethod = "cluster_fast_greedy",
                            weightDeg =  FALSE, normDeg = FALSE,
                            centrLCC = TRUE,
                            hubPar = c("degree", "eigenvector", "betweeness",
                                       "closeness"))
  
  assign(paste0("network_", env), micro_net)
  assign(paste0("properties_", env), net_prop)
  
  edges <- dplyr::select(micro_net$edgelist1, v1, v2)
  edges$Weight <- micro_net$edgelist1$adja
  
  write.csv(edges, file = paste0(working_path, env,"_edges.tsv"), 
            row.names = FALSE)
}



# Initialize an empty list to store the non-null values
properties_list <- list()

# Loop through environments
average_path <- net_prop$globalPropsLCC$avPath1
nodes <- net_prop$globalPropsLCC$lccSize1
net



for (env in environments) {
  prop_env <- get(paste0("properties_", env))
  prop_value <- cbind(prop_value, prop_env[["globalPropsLCC"]])
  
  column_name <- paste0("net_prop$", env)
  properties_list[[column_name]] <- prop_value
  # Write non-null values to CSV file for each environment
}
write.csv(prop_value, paste0(working_path,"net_properties.csv"))
View(properties_list)

# 
# 
# for (env in environments) {
#   prop_env <- get(paste0("properties_", env))
#   prop <- cbind(prop, prop_env[["globalPropsLCC"]])
#   prop[,env] <-  do.call(rbind, prop_env[["globalPropsLCC"]])
#   
# }
# 
# rownames(df) <- row.names(do.call(rbind, prop_env[["globalPropsLCC"]]))
#     