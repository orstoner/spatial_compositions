
# Set up experiment 3 -----------------------------------------------------

# Set the random number generator seed.
set.seed(546237223)

# Number of observed locations.
n_obs_3 <- 1000

# Define the parameter grid.
experiment_3_specs <- expand_grid(
  rep = 1:5,         
  n_basis = 400,
  n_obs = n_obs_3,
  N=c(100,1000,10000,100000),
  chain = c(1,2)
)%>%
  mutate(initial_seed=sample(1:1000000,n()),
         chain_seed=sample(1:1000000,n()))

# Number of chains.
n_chains_3 <- max(experiment_3_specs$chain)

# Number of configurations.
n_specs_3 <- nrow(experiment_3_specs)/n_chains_3

# Unique rep and n_nobs combinations.
observed_rows_3_specs <- experiment_3_specs %>%
  distinct(rep, n_obs) %>%
  arrange(rep, n_obs) %>%
  mutate(dataset_index = row_number())

experiment_3_specs <- experiment_3_specs %>%
  left_join(observed_rows_3_specs, by = c("rep", "n_obs"))

# Create the data sets with
# artificial missing values.
observed_rows_3 <- list()
tree_data_with_missing_3 <- list()

for(i in 1:nrow(observed_rows_3_specs)) {
  n_obs_i <- observed_rows_3_specs$n_obs[i]
  
  # Randomly choose locations to remove 
  # values for all species.
  observed_rows_3[[i]] <- sample(seq_len(n_full), n_obs_i)
  
  tree_data_with_missing_3[[i]] <- tree_model_data%>%
    select(X_coord,Y_coord,larch:sycamore,total)
  tree_data_with_missing_3[[i]][setdiff(seq_len(n_full), observed_rows_3[[i]]),
                                c("larch", "oak", "sitka_spruce", "sycamore")] <- NA
}

# Run experiment 3 --------------------------------------------------------

# Function to run the model for a given configuration.
experiment_3_fun <- function(X) {
  library(dplyr)
  library(magrittr)
  library(nimble)
  
  # Load spatial composition model function.
  source("spatial_composition_model.R")
  
  dataset_index <- experiment_3_specs$dataset_index[X]
  working_data <- tree_data_with_missing_3[[dataset_index]]
  
  # Create the artificial counts.
  working_data%<>%mutate(across(larch:sycamore,
                                ~round(.*experiment_3_specs$N[X])))
  
  # Create total counts for modelling.
  working_totals <- pmax(
    rowSums(select(working_data, larch:sycamore), na.rm = TRUE),
    experiment_3_specs$N[X]
  )
  
  # Run the model.
  run_time <- system.time({
    outputs <- spatial_composition_model(1,
                                         initial_seed = experiment_3_specs$initial_seed[X],
                                         chain_seed = experiment_3_specs$chain_seed[X],
                                         counts_full = select(working_data,larch:sycamore),
                                         locations_full = select(working_data,X_coord,Y_coord),
                                         totals_full = working_totals,
                                         n_types = n_types,
                                         n_basis = experiment_3_specs$n_basis[X],
                                         niter = 100000,
                                         nburn = 60000,
                                         nthin =   40,
                                         define_smooth_using_train = TRUE)
  })
  
  return(list(run_time=run_time,outputs=outputs))
  
}

# Number of runs.
n_runs_3 <- nrow(experiment_3_specs)

# Set up a parallel cluster.
this_cluster <- makeCluster(14)

clusterExport(this_cluster, c("experiment_3_specs", "tree_data_with_missing_3",
                              "spatial_composition_model", "n_types"))

# Run the models.
total_run_time_3 <- system.time({
  experiment_3_outputs <- parLapply(cl = this_cluster, 
                                    X = 1:n_runs_3,
                                    fun = experiment_3_fun)
})

# Add run times to the configuration data frame.
experiment_3_specs$run_time <- do.call("c",lapply(experiment_3_outputs,function(u)u$run_time[3]))

# Assess convergence ------------------------------------------------------

mcmc_lists_3 <- list()
gelman_psrf_3 <- list()
for(i in 1:n_specs_3){
  mcmc_lists_3[[i]] <- as.mcmc.list(lapply(experiment_3_outputs[(n_chains_3*(i-1)+1):(n_chains_3*i)],
                                           function(u)u$outputs$samples))
}

registerDoParallel(cores=14)
gelman_psrf_3 <- foreach(i=1:n_specs_3)%dopar%{
  sapply(1:ncol(mcmc_lists_3[[i]][[1]]),
         function(u)gelman.diag(mcmc_lists_3[[i]][,u],
                                autoburnin = FALSE,
                                multivariate = FALSE,
                                transform = TRUE)$psrf[,1])
}

ess_3 <- foreach(i = 1:n_specs_3, .packages = "coda") %dopar% {
  effectiveSize(mcmc_lists_3[[i]])
}

experiment_3_specs$mean_ess <- rep(
  do.call("c", lapply(ess_3, mean, na.rm = TRUE)),
  each = n_chains_3
)

mean(filter(experiment_3_specs,N==100)$mean_ess)
mean(filter(experiment_3_specs,N==10^5)$mean_ess)

experiment_3_specs$mean_psrf <- rep(do.call("c",lapply(gelman_psrf_3,mean,na.rm=T)),
                                    each=n_chains_3)
experiment_3_specs$p_1.05 <- rep(do.call("c",lapply(gelman_psrf_3,function(u)mean(u<=1.05,na.rm=T))),
                                    each=n_chains_3)
experiment_3_specs$p_1.2 <- rep(do.call("c",lapply(gelman_psrf_3,function(u)mean(u<=1.2,na.rm=T))),
                                    each=n_chains_3)

# Process model outputs ---------------------------------------------------

# Number of MCMC samples per run.
n_samples_3 <- nrow(mcmc_lists_3[[1]][[1]])*n_chains_3

# Posterior predictive samples for all locations.
x_pred_samples_3 <- foreach(i=1:n_specs_3)%dopar%{
  combined_samples <- do.call("rbind",mcmc_lists_3[[i]])
  
  # Coefficient samples.
  b_samples <- subset_b(combined_samples)%>%unlist()%>%
    array(dim=c(n_samples_3,experiment_3_specs$n_basis[n_chains_3*i],n_types))
  
  Z_full <- experiment_3_outputs[[n_chains_3*i]]$outputs$Z_full
  
  mu_samples <- array(NA,dim=c(n_samples_3,n_full,n_types))
  for(d in 1:n_types){
    mu_samples[,,d] <- expit(b_samples[,,d]%*%t(Z_full))
  }
  
  dataset_index <- experiment_3_specs$dataset_index[n_chains_3 * i]
  train_rows <- observed_rows_3[[dataset_index]]
  
  ## Compute variance parameter.
  mu_mean_train <- apply(mu_samples[,train_rows,],c(1,3),mean)
  
  # Polynomial coefficients.
  psi_samples <- subset_psi(combined_samples)%>%unlist()%>%
    array(dim=c(n_samples_3,4,n_types))
  
  phi_samples <- array(NA, dim = c(n_samples_3, n_full, n_types))
  
  for(d in 1:n_types){
    phi_samples[,,d] <- exp(psi_samples[,1, d] +
      psi_samples[,2, d]*(mu_samples[,, d]-mu_mean_train[,d]) +
      psi_samples[,3, d]*(mu_samples[,, d]-mu_mean_train[,d])^2 +
      psi_samples[,4, d]*(mu_samples[,, d]-mu_mean_train[,d])^3)
  }
  
  # Posterior predictive simulation 
  # from the Beta-Binomials.
  x_pred_samples <- array(NA, dim = c(n_samples_3, n_full, n_types))
  
  N_i <- experiment_3_specs$N[n_chains_3 * i]
  
  remainder <- matrix(N_i, nrow = n_samples_3, ncol = n_full)
  
  for(d in seq_len(n_types)) {
    
    p_d <- rbeta(
      n_samples_3 * n_full,
      shape1 = as.vector(mu_samples[,,d] * phi_samples[,,d]),
      shape2 = as.vector((1 - mu_samples[,,d]) * phi_samples[,,d])
    )
    
    x_d <- rbinom(
      n_samples_3 * n_full,
      size = as.vector(remainder),
      prob = p_d
    )
    
    x_d <- matrix(x_d, nrow = n_samples_3, ncol = n_full)
    
    x_pred_samples[,,d] <- x_d
    
    remainder <- remainder - x_d
  }
  
  x_pred_samples/N_i
}

# Compute posterior predictive summaries.
pred_summaries_3 <-  foreach(i=1:n_specs_3)%dopar%{
  
  # Reload the original data.
  tree_proportions <- readRDS("Data/tree_proportions.rds")%>%
    select(X_coord,Y_coord,larch,oak,sitka_spruce,sycamore,total)
  
  tree_proportions$observed=FALSE
  tree_proportions$observed[observed_rows_3[[experiment_3_specs$dataset_index[n_chains_3*i]]]] <- TRUE
  
  # Original complete values as a matrix.
  actual_mat <- tree_proportions %>%
    select(larch:sycamore) %>%
    as.matrix()
  
  # CRPS matrix: n_full x n_types
  crps_mat <- matrix(NA,nrow=n_full,ncol=n_types)
  
  for (j in 1:n_full) {
    for (d in 1:n_types) {
      crps_mat[j, d] <- crps_from_samples(
        samples = x_pred_samples_3[[i]][, j, d],
        y = actual_mat[j, d]
      )
    }
  }
  
  # Put all summaries together in long format.
  tree_proportions_long <- tree_proportions%>%pivot_longer(larch:sycamore,names_to="species",
                                                   values_to="actual")%>%
    mutate(lower=as.numeric(t(apply(x_pred_samples_3[[i]],c(2,3),quantile,0.025))),
           median=as.numeric(t(apply(x_pred_samples_3[[i]],c(2,3),median))),
           mean=as.numeric(t(apply(x_pred_samples_3[[i]],c(2,3),mean))),
           upper=as.numeric(t(apply(x_pred_samples_3[[i]],c(2,3),quantile,0.975))),
           crps   = as.numeric(t(crps_mat)),
           n_basis=experiment_3_specs$n_basis[n_chains_3*i],
           n_obs=experiment_3_specs$n_obs[n_chains_3*i],
           rep=experiment_3_specs$rep[n_chains_3*i],
           N=experiment_3_specs$N[n_chains_3*i])%>%
    group_by(species)%>%
    mutate(naive=mean(actual[observed==TRUE]))
}

pred_summaries_3%<>%do.call("rbind",.)

# Assess prediction performance -------------------------------------------

# Compute performance summaries.
pred_perf_3 <- pred_summaries_3%>%
  group_by(N,observed,rep)%>%
  summarise(mae=mean(abs(actual-mean)),
            rmse=sqrt(mean((actual-mean)^2)),
            mean_crps=mean(crps),
            coverage=mean(actual>=lower&actual<=upper),
            xi=1-mean((actual-mean)^2)/mean((actual-naive)^2))

pred_perf_3 


# Plot of mean absolute error and mean CRPS
# by artificial total count N.
N_plot <- ggplot(pred_perf_3%>%
                         select(N,rep,observed,mean_crps,mae,rmse)%>%
                         rename(`Mean CRPS`=mean_crps,
                                `Mean absolute error`=mae,
                                `Root mean square error`=rmse)%>%
                         pivot_longer(cols=c("Mean CRPS","Mean absolute error","Root mean square error"),
                                      names_to = "metric")%>%
                         mutate(observed=case_match(observed,
                                                    FALSE ~ "Unobserved locations",
                                                    TRUE ~ "Observed locations"),
                                metric=factor(metric,levels=c("Mean absolute error","Root mean square error",
                                                              "Mean CRPS"))),
                       aes(x=N,y=value,colour=observed))+
  facet_wrap(~metric,scales="free")+
  geom_smooth(linewidth=0.5,method="gam",formula=y ~s(x,k=4))+
  #geom_jitter(height=0)+
  geom_point()+
  scale_colour_viridis_d(begin = 0.1,end=0.9)+
  theme(panel.grid = element_blank())+
  labs(x="Assumed total count (N)",y=NULL,
       colour=NULL)+
  scale_shape_manual(values=15:18)+
  scale_x_log10()
ggsave(N_plot,filename="Plots/N.pdf",width=8,height=2.25)


ggplot(pred_perf_3%>%
         left_join(select(experiment_3_specs,rep,N,p_1.05)%>%distinct()),
       aes(x=N,y=rmse,colour=p_1.05))+
  geom_point()+
  scale_x_log10()

ggplot(pred_perf_3,
       aes(x=N,y=rmse,colour=observed))+
  geom_point()+
  scale_x_log10()

ggplot(pred_perf_3,
       aes(x=N,y=mae,colour=observed))+
  geom_point()+
  scale_x_log10()


ggplot(experiment_3_specs,
       aes(x=N,y=run_time))+
  geom_point()+
  scale_x_log10()
# Check zero proportions --------------------------------------------------

# Context labels for experiment 3.
location_context_labels <- c(
  "Observed locations",
  "Unobserved locations"
)

# Compute posterior predictive proportion of locations with zeros
# by species, context, N, replicate, and MCMC sample.
zero_props_3 <- foreach(i = seq_len(n_specs_3), .combine = bind_rows) %dopar% {
  
  dataset_index <- experiment_3_specs$dataset_index[n_chains_3 * i]
  
  samples_i <- x_pred_samples_3[[i]]
  # dim: n_samples_3 x n_full x n_types
  
  observed_i <- rep(FALSE, n_full)
  observed_i[observed_rows_3[[dataset_index]]] <- TRUE
  
  context_i <- ifelse(
    observed_i,
    "Observed locations",
    "Unobserved locations"
  )
  
  out <- expand.grid(
    sample = seq_len(dim(samples_i)[1]),
    species = tree_types,
    context = location_context_labels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    mutate(
      species_index = match(species, tree_types),
      zero_prop = purrr::pmap_dbl(
        list(sample, species_index, context),
        function(s, d, ctx) {
          
          rows_ctx <- which(context_i == ctx)
          
          if(length(rows_ctx) == 0) {
            return(NA_real_)
          }
          
          mean(samples_i[s, rows_ctx, d] == 0, na.rm = TRUE)
        }
      ),
      n_basis = experiment_3_specs$n_basis[n_chains_3 * i],
      n_obs = experiment_3_specs$n_obs[n_chains_3 * i],
      N = experiment_3_specs$N[n_chains_3 * i],
      rep = experiment_3_specs$rep[n_chains_3 * i],
      dataset_index = dataset_index
    ) %>%
    select(
      sample,
      species,
      context,
      zero_prop,
      n_basis,
      n_obs,
      N,
      rep,
      dataset_index
    )
  
  out
}

# Reload the original data.
truth_full <- readRDS("Data/tree_proportions.rds") %>%
  select(all_of(tree_types)) %>%
  as.matrix()

# Compute true proportion of locations with zeros
# by species and context.
zero_props_truth_3 <- foreach(i = seq_len(n_specs_3), .combine = bind_rows) %dopar% {
  
  dataset_index <- experiment_3_specs$dataset_index[n_chains_3 * i]
  
  observed_i <- rep(FALSE, n_full)
  observed_i[observed_rows_3[[dataset_index]]] <- TRUE
  
  context_i <- ifelse(
    observed_i,
    "Observed locations",
    "Unobserved locations"
  )
  
  out <- expand.grid(
    species = tree_types,
    context = location_context_labels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    mutate(
      species_index = match(species, tree_types),
      zero_prop = purrr::map2_dbl(
        species_index,
        context,
        function(d, ctx) {
          
          rows_ctx <- which(context_i == ctx)
          
          if(length(rows_ctx) == 0) {
            return(NA_real_)
          }
          
          mean(truth_full[rows_ctx, d] == 0, na.rm = TRUE)
        }
      ),
      n_basis = experiment_3_specs$n_basis[n_chains_3 * i],
      n_obs = experiment_3_specs$n_obs[n_chains_3 * i],
      N = experiment_3_specs$N[n_chains_3 * i],
      rep = experiment_3_specs$rep[n_chains_3 * i],
      dataset_index = dataset_index
    ) %>%
    select(
      species,
      context,
      zero_prop,
      n_basis,
      n_obs,
      N,
      rep,
      dataset_index
    )
  
  out
}

# Summarise posterior predictive zero proportion distributions.
zero_props_summary_3 <- zero_props_3 %>%
  group_by(species, context, N, rep, n_basis, n_obs, dataset_index) %>%
  summarise(
    lower = quantile(zero_prop, 0.025, na.rm = TRUE),
    mean = mean(zero_prop, na.rm = TRUE),
    upper = quantile(zero_prop, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    zero_props_truth_3,
    by = c(
      "species",
      "context",
      "N",
      "rep",
      "n_basis",
      "n_obs",
      "dataset_index"
    ),
    suffix = c("_pred", "_truth")
  ) %>%
  mutate(
    context = factor(
      context,
      levels = location_context_labels
    ),
    species = case_when(
      species == "larch" ~ "Larch",
      species == "oak" ~ "Oak",
      species == "sitka_spruce" ~ "Sitka spruce",
      species == "sycamore" ~ "Sycamore",
      TRUE ~ species
    )
  )
# Plot of posterior predictive zero proportion
# compared to the original data.
N_zero_plot <- ggplot(
  zero_props_summary_3,
  aes(x = mean, y = zero_prop, shape = species, colour = species)
) +
  facet_grid(context ~ N, labeller=labeller(
    N = as_labeller(function(x) paste("N", x, sep = " = ")),
    context = label_value
  )) +
  geom_abline() +
  geom_errorbar(aes(xmin = lower, xmax = upper)) +
  geom_point() +
  theme(panel.grid = element_blank()) +
  labs(
    x = "Predicted zero proportion",
    y = "Observed zero proportion",
    colour = NULL,
    shape = NULL
  ) +
  scale_shape_manual(values = c(15:18)) +
  scale_colour_viridis_d(begin = 0.1, end = 0.9)
ggsave(N_zero_plot,filename="Plots/N_zero_prop.pdf",width=8,height=3.5)

# Zero proportion 95% PI coverage.
zero_props_summary_3%>%
  group_by(N,context)%>%
  summarise(mean(zero_prop>=lower&zero_prop<=upper))

# Zero proportion MAE.
zero_props_summary_3%>%
  group_by(N,context)%>%
  summarise(mean(abs(zero_prop-mean)))

# Zero proportion bias.
zero_props_summary_3%>%
  group_by(N,context)%>%
  summarise(mean((zero_prop-mean)))

select(experiment_3_specs,rep,N,run_time)%>%distinct()

# Correlation between log-10 N and run-time.
cor(log((select(experiment_3_specs,rep,N,run_time)%>%distinct())$N,10),
    (select(experiment_3_specs,rep,N,run_time)%>%distinct())$run_time)

# Save outputs ------------------------------------------------------------

# Save pretty much everything.
save(experiment_3_specs,observed_rows_3,
     n_specs_3,n_chains_3,n_samples_3,n_obs_3,
     total_run_time_3,experiment_3_outputs,
     tree_data_with_missing_3,
     zero_props_3,zero_props_summary_3,
     zero_props_truth_3,
     x_pred_samples_3,pred_summaries_3,
     pred_perf_3,file="Saved runs/experiment_3_12_06_36.RData")

# Save only smaller objects required to inspect
# and plot the results.
save(experiment_3_specs,observed_rows_3,
     n_specs_3,n_chains_3,n_samples_3,n_obs_3,
     total_run_time_3,
     tree_data_with_missing_3,
     zero_props_3,zero_props_summary_3,
     zero_props_truth_3,
     pred_summaries_3,
     pred_perf_3,file="Saved runs/experiment_3_lite_12_06_36.RData")
