
# Load required packages --------------------------------------------------

library(dplyr)
library(tidyr)
library(nimble)
library(nimbleHMC)
library(coda)
library(ggplot2)
library(mgcv)
library(parallel)
library(magrittr)
library(doParallel)
library(glmmTMB)
library(openxlsx)
library(reshape2)
library(compositions)

# Path to the tree proportions data.
tree_proportions_path <- "Data/simulated_tree_proportions.rds"

# Load spatial composition model function.
source("spatial_composition_model.R")

# Prepare data for modelling ----------------------------------------------

# Load the data.
tree_proportions <- readRDS(tree_proportions_path)

# Heatmap of the original data.
png("Plots/heatmap_originals.png", height = 4,width=8,units = "in",res=600)
plot(ggplot(tree_proportions%>%
              pivot_longer(cols =larch:sycamore,names_to="type",values_to = "value")%>%
              mutate(type=case_match(type,"larch"~"Larch",
                                     "oak"~"Oak","scots_pine"~"Scots pine",
                                     "shadow"~"Shadow","silver_birch"~"Silver birch",
                                     "sitka_spruce"~"Sitka spruce","sweet_chestnut"~"Sweet chestnut",
                                     "sycamore"~"Sycamore")),
            aes(x=X_coord,y=Y_coord,fill=value/total))+
       geom_tile()+
       facet_wrap(~type,nrow=2)+
       scale_fill_viridis_c(name="Proportion",breaks=seq(0,1,by=0.2)) +
       scale_colour_viridis_c(name="Proportion",breaks=seq(0,1,by=0.2)) +
       theme(axis.text.x=element_text(angle=60, hjust=1)) +
       coord_fixed()+
       labs(x = "", y = "", fill = 'Count') +
       theme( axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.grid = element_blank(),
              legend.position = "bottom",
              legend.key.width = unit(0.5, "in"),
              legend.title = element_text(vjust=0.8)))
dev.off()

# Four species names.
tree_types <- c("larch",
                "oak",
                "sitka_spruce",
                "sycamore")

# Reduce to four species.
tree_proportions <- tree_proportions%>%
  select(X_coord,Y_coord,tree_types,total)

# Number of species.
n_types <- length(tree_types)

# Center and scale x and y coordinates.
tree_model_data <- tree_proportions%>% 
  select(X_coord,
         Y_coord,
         all_of(tree_types),
         total) 

# Full number of locations.
n_full <- nrow(tree_model_data)

# Center.
tree_model_data$X_coord <- tree_model_data$X_coord-mean(tree_model_data$X_coord)
tree_model_data$Y_coord <- tree_model_data$Y_coord-mean(tree_model_data$Y_coord)

# Scale by maximum pairwise distance.
max_dist <- max(dist(as.matrix(tree_model_data[,c("X_coord","Y_coord")])))

tree_model_data$X_coord <- tree_model_data$X_coord/max_dist
tree_model_data$Y_coord <- tree_model_data$Y_coord/max_dist



