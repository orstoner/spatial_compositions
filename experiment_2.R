
# Set up experiment 2 -----------------------------------------------------

# Set the random number generator seed.
set.seed(982348432)

# Number of observed locations.
n_obs_2 <- c(100,200,400)

# Define the parameter grid.
experiment_2_specs <- expand_grid(
  rep = 1:5,         
  n_basis = c(50,100,200,400,600),
  n_obs = n_obs_2,
  N=100,
  chain = c(1,2)
)%>%
  mutate(initial_seed=sample(1:1000000,n()),
         chain_seed=sample(1:1000000,n()))

# Number of chains.
n_chains_2 <- max(experiment_2_specs$chain)

# Number of configurations.
n_specs_2 <- nrow(experiment_2_specs)/n_chains_2

# Unique rep and n_nobs combinations.
observed_rows_2_specs <- experiment_2_specs %>%
  distinct(rep, n_obs) %>%
  arrange(rep, n_obs) %>%
  mutate(dataset_index = row_number())

experiment_2_specs <- experiment_2_specs %>%
  left_join(observed_rows_2_specs, by = c("rep", "n_obs"))

# Create the data sets with
# artificial missing values.
observed_rows_2 <- list()
tree_data_with_missing_2 <- list()

for(i in 1:nrow(observed_rows_2_specs)) {
  n_obs_i <- observed_rows_2_specs$n_obs[i]
  
  # Randomly choose locations to remove 
  # values for all species.
  observed_rows_2[[i]] <- sample(seq_len(n_full), n_obs_i)
  
  tree_data_with_missing_2[[i]] <- tree_model_data
  tree_data_with_missing_2[[i]][setdiff(seq_len(n_full), observed_rows_2[[i]]),
                                c("larch", "oak", "sitka_spruce", "sycamore")] <- NA
}

# Run experiment 2 --------------------------------------------------------

# Function to run the model for a given configuration.
experiment_2_fun <- function(X) {
  library(dplyr)
  library(magrittr)
  library(nimble)
  
  # Load spatial composition model function.
  source("spatial_composition_model.R")
  
  dataset_index <- experiment_2_specs$dataset_index[X]
  working_data <- tree_data_with_missing_2[[dataset_index]]
  
  # Create the artificial counts.
  working_data%<>%mutate(across(larch:sycamore,
                                ~round(.*experiment_2_specs$N[X])))
  
  # Create total counts for modelling.
  working_totals <- pmax(
    rowSums(select(working_data, larch:sycamore), na.rm = TRUE),
    experiment_2_specs$N[X]
  )
  
  # Run the model.
  run_time <- system.time({
    outputs <- spatial_composition_model(1,
                                         initial_seed = experiment_2_specs$initial_seed[X],
                                         chain_seed = experiment_2_specs$chain_seed[X],
                                         counts_full = select(working_data,larch:sycamore),
                                         locations_full = select(working_data,X_coord,Y_coord),
                                         totals_full = working_totals,
                                         n_types = n_types,
                                         n_basis = experiment_2_specs$n_basis[X],
                                         niter = 100000,
                                         nburn = 60000,
                                         nthin =   40,
                                         define_smooth_using_train = FALSE)
  })
  
  return(list(run_time=run_time,outputs=outputs))
  
}

# Number of runs.
n_runs_2 <- nrow(experiment_2_specs)

# Set up a parallel cluster.
this_cluster <- makeCluster(14)

clusterExport(this_cluster, c("experiment_2_specs", "tree_data_with_missing_2",
                              "spatial_composition_model", "n_types"))

# Run the models.
total_run_time_2 <- system.time({
  experiment_2_outputs <- parLapply(cl = this_cluster, 
                                    X = 1:n_runs_2,
                                    fun = experiment_2_fun)
})

# Add run times to the configuration data frame.
experiment_2_specs$run_time <- do.call("c",lapply(experiment_2_outputs,function(u)u$run_time[3]))

# Assess convergence ------------------------------------------------------

mcmc_lists_2 <- list()
gelman_psrf_2 <- list()
for(i in 1:n_specs_2){
  mcmc_lists_2[[i]] <- as.mcmc.list(lapply(experiment_2_outputs[(n_chains_2*(i-1)+1):(n_chains_2*i)],
                                           function(u)u$outputs$samples))
}

registerDoParallel(cores=14)
gelman_psrf_2 <- foreach(i=1:n_specs_2)%dopar%{
  sapply(1:ncol(mcmc_lists_2[[i]][[1]]),
         function(u)gelman.diag(mcmc_lists_2[[i]][,u],
                                autoburnin = FALSE,
                                multivariate = FALSE,
                                transform = TRUE)$psrf[,1])
}

experiment_2_specs$mean_psrf <- rep(do.call("c",lapply(gelman_psrf_2,mean,na.rm=T)),
                                    each=n_chains_2)
experiment_2_specs$p_1.05 <- rep(do.call("c",lapply(gelman_psrf_2,function(u)mean(u<=1.05,na.rm=T))),
                                    each=n_chains_2)
experiment_2_specs$p_1.2 <- rep(do.call("c",lapply(gelman_psrf_2,function(u)mean(u<=1.2,na.rm=T))),
                                    each=n_chains_2)

# Process model outputs ---------------------------------------------------

# Number of MCMC samples per run.
n_samples_2 <- nrow(mcmc_lists_2[[1]][[1]])*n_chains_2

# Posterior predictive samples for all locations.
x_pred_samples_2 <- foreach(i=1:n_specs_2)%dopar%{
  combined_samples <- do.call("rbind",mcmc_lists_2[[i]])
  
  # Coefficient samples.
  b_samples <- subset_b(combined_samples)%>%unlist()%>%
    array(dim=c(n_samples_2,experiment_2_specs$n_basis[n_chains_2*i],n_types))
  
  Z_full <- experiment_2_outputs[[n_chains_2*i]]$outputs$Z_full
  
  mu_samples <- array(NA,dim=c(n_samples_2,n_full,n_types))
  for(d in 1:n_types){
    mu_samples[,,d] <- expit(b_samples[,,d]%*%t(Z_full))
  }
  
  dataset_index <- experiment_2_specs$dataset_index[n_chains_2 * i]
  train_rows <- observed_rows_2[[dataset_index]]
  
  ## Compute variance parameter.
  mu_mean_train <- apply(mu_samples[,train_rows,],c(1,3),mean)
  
  # Polynomial coefficients.
  psi_samples <- subset_psi(combined_samples)%>%unlist()%>%
    array(dim=c(n_samples_2,4,n_types))
  
  phi_samples <- array(NA, dim = c(n_samples_2, n_full, n_types))
  
  for(d in 1:n_types){
    phi_samples[,,d] <- exp(psi_samples[,1, d] +
                              psi_samples[,2, d]*(mu_samples[,, d]-mu_mean_train[,d]) +
                              psi_samples[,3, d]*(mu_samples[,, d]-mu_mean_train[,d])^2 +
                              psi_samples[,4, d]*(mu_samples[,, d]-mu_mean_train[,d])^3)
  }
  
  # Posterior predictive simulation 
  # from the Beta-Binomials.
  x_pred_samples <- array(NA, dim = c(n_samples_2, n_full, n_types))
  
  N_i <- experiment_2_specs$N[n_chains_2 * i]
  
  
  remainder <- matrix(N_i, nrow = n_samples_2, ncol = n_full)
  
  for(d in seq_len(n_types)) {
    
    p_d <- rbeta(
      n_samples_2 * n_full,
      shape1 = as.vector(mu_samples[,,d] * phi_samples[,,d]),
      shape2 = as.vector((1 - mu_samples[,,d]) * phi_samples[,,d])
    )
    
    x_d <- rbinom(
      n_samples_2 * n_full,
      size = as.vector(remainder),
      prob = p_d
    )
    
    x_d <- matrix(x_d, nrow = n_samples_2, ncol = n_full)
    
    x_pred_samples[,,d] <- x_d
    
    remainder <- remainder - x_d
  }
  
  x_pred_samples/N_i
}

# Compute posterior predictive summaries.
pred_summaries_2 <-  foreach(i=1:n_specs_2)%dopar%{
  
  # Reload the original data.
  tree_proportions <- readRDS(tree_proportions_path)%>%
    select(X_coord,Y_coord,larch,oak,sitka_spruce,sycamore)
  
  tree_proportions$observed=FALSE
  tree_proportions$observed[observed_rows_2[[experiment_2_specs$dataset_index[n_chains_2*i]]]] <- TRUE
  
  # Original complete values as a matrix.
  actual_mat <- tree_proportions %>%
    select(larch:sycamore) %>%
    as.matrix()
  
  # CRPS matrix: n_full x n_types
  crps_mat <- matrix(NA,nrow=n_full,ncol=n_types)
  
  for (j in 1:n_full) {
    for (d in 1:n_types) {
      crps_mat[j, d] <- crps_from_samples(
        samples = x_pred_samples_2[[i]][, j, d],
        y = actual_mat[j, d]
      )
    }
  }
  
  # Put all summaries together in long format.
  tree_proportions_long <- tree_proportions%>%pivot_longer(larch:sycamore,names_to="species",
                                                   values_to="actual")%>%
    mutate(lower=as.numeric(t(apply(x_pred_samples_2[[i]],c(2,3),quantile,0.025))),
           median=as.numeric(t(apply(x_pred_samples_2[[i]],c(2,3),median))),
           mean=as.numeric(t(apply(x_pred_samples_2[[i]],c(2,3),mean))),
           upper=as.numeric(t(apply(x_pred_samples_2[[i]],c(2,3),quantile,0.975))),
           crps   = as.numeric(t(crps_mat)),
           n_basis=experiment_2_specs$n_basis[n_chains_2*i],
           n_obs=experiment_2_specs$n_obs[n_chains_2*i],
           rep=experiment_2_specs$rep[n_chains_2*i])%>%
    group_by(species)%>%
    mutate(naive=mean(actual[observed==TRUE]))
}

pred_summaries_2%<>%do.call("rbind",.)

# Assess prediction performance -------------------------------------------

# Compute performance summaries.
pred_perf_2 <- pred_summaries_2%>%
  group_by(n_obs,n_basis,observed,rep)%>%
  summarise(mae=mean(abs(actual-mean)),
            mean_crps=mean(crps),
            coverage=mean(actual>=lower&actual<=upper),
            xi=1-mean((actual-mean)^2)/mean((actual-naive)^2))

pred_perf_2 

# Plot of mean absolute error and mean CRPS
# by number of basis functions.
n_basis_plot <- ggplot(pred_perf_2%>%
         select(n_basis,rep,observed,n_obs,mean_crps,mae)%>%
         rename(`Mean CRPS`=mean_crps,
                `Mean absolute error`=mae)%>%
         pivot_longer(cols=c("Mean CRPS","Mean absolute error"),
                      names_to = "metric")%>%
         mutate(observed=case_match(observed,
                                    FALSE ~ "Unobserved locations",
                                    TRUE ~ "Observed locations")),
       aes(x=n_basis,y=value,colour=observed))+
  facet_wrap(~metric,scales="free")+
  geom_smooth(linewidth=0.5,method="gam",formula=y ~s(x,k=5))+
  geom_jitter(aes(shape=factor(n_obs,levels=rev(n_obs_2))),height=0)+
  scale_colour_viridis_d(begin = 0.1,end=0.9)+
  theme(panel.grid = element_blank())+
  labs(x="Number of basis functions",y=NULL,
       shape="Number of observed\nlocations",
       colour=NULL)+
  scale_shape_manual(values=15:18)
ggsave(n_basis_plot,filename="Plots/n_basis.pdf",width=8,height=3)

# Plot of run time
# by number of basis functions.
run_time_plot <- ggplot(experiment_2_specs,
       aes(x=n_basis,y=run_time/60,
           shape=factor(n_obs,levels=rev(n_obs_2)),
           colour=factor(n_obs,levels=rev(n_obs_2))))+
  geom_smooth(linewidth=0.5,method="gam",formula=y~s(x,k=5))+
  geom_jitter(height=0)+
  theme(panel.grid = element_blank())+
  labs(x="Number of basis functions",y="Computation time (minutes)",
       shape="Number of\nobserved\nlocations",
       colour="Number of\nobserved\nlocations")+
  scale_shape_manual(values=15:19)+
  scale_colour_viridis_d(begin=0.1,end=0.9,direction=-1)
ggsave(run_time_plot,filename="Plots/run_time.pdf",width=4,height=2.5)

# Save outputs ------------------------------------------------------------

# Save pretty much everything.
save(experiment_2_specs,observed_rows_2,
     n_specs_2,n_chains_2,n_samples_2,n_obs_2,
     total_run_time_2,experiment_2_outputs,
     tree_data_with_missing_2,
     x_pred_samples_2,pred_summaries_2,
     pred_perf_2,file="Saved runs/experiment_2_14_06_26.RData")

# Save only smaller objects required to inspect
# and plot the results.
save(experiment_2_specs,observed_rows_2,
     n_specs_2,n_chains_2,n_samples_2,n_obs_2,
     total_run_time_2,
     tree_data_with_missing_2,
     pred_summaries_2,
     pred_perf_2,file="Saved runs/experiment_2_lite_14_06_26.RData")

