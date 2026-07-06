
# Create square grid ------------------------------------------------------

# Length of each side.
n_side <- 46

sim_grid <- expand.grid(
  X_coord = seq(0, 1, length.out = n_side),
  Y_coord = seq(0, 1, length.out = n_side)
)

# Center.
sim_grid$X_coord <- sim_grid$X_coord-mean(sim_grid$X_coord)
sim_grid$Y_coord <- sim_grid$Y_coord-mean(sim_grid$Y_coord)

# Scale by maximum pairwise distance.
max_dist <- max(dist(as.matrix(sim_grid[,c("X_coord","Y_coord")])))

sim_grid$X_coord <- sim_grid$X_coord/max_dist
sim_grid$Y_coord <- sim_grid$Y_coord/max_dist

n_locations <- nrow(sim_grid)
n_components <- 10  

# Simulate proportions -------------------------------------------------------

sim_data_matrix <- matrix(NA, nrow = n_locations, ncol = n_components)

# Set up spatial basis functions.
sim_basis <- smoothCon(
  s(X_coord, Y_coord, bs = 'gp', m = c(-2, 0.1, 1), k = n_locations/2),
  data = sim_grid,
  absorb.cons = TRUE
)[[1]]

sim_X <- PredictMat(sim_basis,sim_grid)

# Simulate for each component.
set.seed(346578238)
for (i in 1:n_components) {
  # Simulate smooth coefficients
  sim_beta <- rnorm(ncol(sim_X), sd = 1)
  
  # Store simulated raw values
  sim_data_matrix[, i] <- rnorm(1,0,0.5) +  # Intercept.
    2*as.vector(scale(sim_X %*% sim_beta)) + # Latent effect.
    rnorm(n_locations,0,0.1) # Noise term.
}

# Convert to proportions using the inverse CLR transformation.

proportions <- t(apply(sim_data_matrix, 1, function(row) {
  inv_clr <- exp(row) / sum(exp(row))
  return(round(inv_clr,2)) # Round to 2 decimal places.
}))

# Output the simulated data -----------------------------------------------

# Format as a data frame.
simulated_data <- data.frame(
  sim_grid,
  proportions
)

# Rename columns to match the original data.
names(simulated_data) <- c('X_coord', 'Y_coord','ash', 'beech', 'larch', 'oak', 'scots_pine', 
                           'shadow', 'silver_birch', 'sitka_spruce', 'sweet_chestnut', 'sycamore')

# Order columns to match the original data.
simulated_data <- simulated_data %>%
  dplyr::select(ash:sycamore,
         ends_with('_coord'))%>%
  mutate(total=rowSums(across(ash:sycamore)))

# Heatmap of the simulated data.
png("Plots/heatmap_simulated.png", height = 4,width=8,units = "in",res=600)
plot(ggplot(simulated_data%>%
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

# Save the data.
saveRDS(simulated_data, "Data/simulated_tree_proportions.rds")
