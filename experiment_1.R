
# Set up experiment 1 -----------------------------------------------------

# Set the random number generator seed.
set.seed(982348432)

# Number of observed locations.
n_obs_1 <- 1000

# Define the parameter grid.
experiment_1_specs <- expand_grid(
  rep = 1:10,         
  n_basis = 400,
  n_obs = n_obs_1,
  N=100,
  include_polynomial=c(TRUE,FALSE),
  chain = c(1,2)
)%>%
  mutate(initial_seed=sample(1:1000000,n()),
         chain_seed=sample(1:1000000,n()))

# Number of chains.
n_chains_1 <- max(experiment_1_specs$chain)

# Number of configurations.
n_specs_1 <- nrow(experiment_1_specs)/n_chains_1

# Unique rep and n_nobs combinations.
observed_rows_1_specs <- experiment_1_specs %>%
  distinct(rep, n_obs) %>%
  arrange(rep, n_obs) %>%
  mutate(dataset_index = row_number())

experiment_1_specs <- experiment_1_specs %>%
  left_join(observed_rows_1_specs, by = c("rep", "n_obs"))

# Create the data sets with
# artificial missing values.
observed_rows_1 <- list()
tree_data_with_missing_1 <- list()

for(i in seq_len(nrow(observed_rows_1_specs))) {
  
  n_obs_i <- observed_rows_1_specs$n_obs[i]
  
  # Locations with at least one species observed.
  observed_rows_1[[i]] <- sample(seq_len(n_full), n_obs_i)
  obs_rows_i <- observed_rows_1[[i]]
  
  # Start from complete data.
  tree_data_with_missing_1[[i]] <- tree_model_data
  
  # Rows outside obs_rows_i have all four species missing.
  n_missing <- rep(length(tree_types), n_full)
  
  # Rows inside obs_rows_i initially have no missing species.
  n_missing[obs_rows_i] <- 0
  
  # Locations with one species missing.
  n_missing_1 <- 200
  
  missing_1_rows <- sample(
    obs_rows_i,
    n_missing_1,
    replace = FALSE
  )
  
  remaining_rows <- setdiff(obs_rows_i, missing_1_rows)
  
  # Locations with two species missing.
  n_missing_2 <- 100
  
  missing_2_rows <- sample(
    remaining_rows,
    n_missing_2,
    replace = FALSE
  )
  
  remaining_rows <- setdiff(obs_rows_i, c(missing_1_rows, missing_2_rows))
  
  # Locations with three species missing.
  n_missing_3 <- min(100, length(remaining_rows))
  
  missing_3_rows <- sample(
    remaining_rows,
    n_missing_3,
    replace = FALSE
  )
  
  n_missing[missing_1_rows] <- 1
  n_missing[missing_2_rows] <- 2
  n_missing[missing_3_rows] <- 3
  
  # Add missing values to random species.
  for(j in seq_len(n_full)) {
    if(n_missing[j] > 0) {
      cols_to_make_missing <- sample(
        tree_types,
        n_missing[j],
        replace = FALSE
      )
      
      tree_data_with_missing_1[[i]][j, cols_to_make_missing] <- NA
    }
  }
  
  # Add missing value identified columns.
  tree_data_with_missing_1[[i]] <- tree_data_with_missing_1[[i]] %>%
    mutate(n_missing = n_missing) %>%
    mutate(
      across(
        all_of(tree_types),
        ~ is.na(.),
        .names = "{.col}_missing"
      )
    )
}


# Run experiment 1 --------------------------------------------------------

# Function to run the model for a given configuration.
experiment_1_fun <- function(X) {
  library(dplyr)
  library(magrittr)
  library(nimble)
  
  # Load spatial composition model function.
  source("spatial_composition_model.R")
  
  dataset_index <- experiment_1_specs$dataset_index[X]
  working_data <- tree_data_with_missing_1[[dataset_index]]
  
  # Create the artificial counts.
  working_data%<>%mutate(across(larch:sycamore,
                                ~round(.*experiment_1_specs$N[X])))
  
  # Create total counts for modelling.
  working_totals <- pmax(
    rowSums(select(working_data, larch:sycamore), na.rm = TRUE),
    experiment_1_specs$N[X]
  )

  # Determine whether to run the model with or without the polynomial.
  include_polynomial <- 
    !"include_polynomial" %in% names(experiment_1_specs) ||
    isTRUE(experiment_1_specs$include_polynomial[X])
  
  # Run the model.
  run_time <- system.time({
    model_fun <- if (include_polynomial) {
      spatial_composition_model
    } else {
      spatial_composition_model_no_polynomial
    }
    
    outputs <- model_fun(
      1,
      initial_seed = experiment_1_specs$initial_seed[X],
      chain_seed   = experiment_1_specs$chain_seed[X],
      counts_full  = working_data %>% select(larch:sycamore),
      locations_full = working_data %>% select(X_coord, Y_coord),
      totals_full            = working_totals,
      n_types      = n_types,
      n_basis      = experiment_1_specs$n_basis[X],
      niter        = 20000,
      nburn        = 12000,
      nthin        = 8,
      define_smooth_using_train = TRUE
    )
    
  })
  
  return(list(run_time=run_time,outputs=outputs))
  
}

# Number of runs.
n_runs_1 <- nrow(experiment_1_specs)

# Set up a parallel cluster.
this_cluster <- makeCluster(10)

clusterExport(this_cluster, c("experiment_1_specs", "tree_data_with_missing_1",
                              "n_types"))

# Run the models.
total_run_time_1 <- system.time({
  experiment_1_outputs <- parLapply(cl = this_cluster, 
                                    X = 1:n_runs_1,
                                    fun = experiment_1_fun)
})

stopCluster(this_cluster)

save(total_run_time_1,experiment_1_outputs,file="Saved runs/experiment_1_outputs.RData")

# Add run times to the configuration data frame.
experiment_1_specs$run_time <- do.call("c",lapply(experiment_1_outputs,function(u)u$run_time[3]))

# Compare mean run time of GDM+ and GDM models.
mean_run_time <- experiment_1_specs%>%group_by(include_polynomial)%>%summarise(mean(run_time))
mean_run_time[2,2]/mean_run_time[1,2]


# Assess convergence ------------------------------------------------------

mcmc_lists_1 <- list()
gelman_psrf_1 <- list()
for(i in 1:n_specs_1){
  mcmc_lists_1[[i]] <- as.mcmc.list(lapply(experiment_1_outputs[(n_chains_1*(i-1)+1):(n_chains_1*i)],
                                           function(u)u$outputs$samples))
}

registerDoParallel(cores=10)
gelman_psrf_1 <- foreach(i=1:n_specs_1)%dopar%{
  sapply(1:ncol(mcmc_lists_1[[i]][[1]]),
         function(u)gelman.diag(mcmc_lists_1[[i]][,u],
                                autoburnin = FALSE,
                                multivariate = FALSE,
                                transform = TRUE)$psrf[,1])
}

experiment_1_specs$mean_psrf <- rep(do.call("c",lapply(gelman_psrf_1,mean,na.rm=T)),
                                    each=n_chains_1)
experiment_1_specs$p_1.05 <- rep(do.call("c",lapply(gelman_psrf_1,function(u)mean(u<=1.05,na.rm=T))),
                                    each=n_chains_1)
experiment_1_specs$p_1.2 <- rep(do.call("c",lapply(gelman_psrf_1,function(u)mean(u<=1.2,na.rm=T))),
                                    each=n_chains_1)


# Process model outputs ---------------------------------------------------

# Number of MCMC samples per run.
n_samples_1 <- nrow(mcmc_lists_1[[1]][[1]])*n_chains_1

# Samples direct from MCMC (training locations only).
x_samples_1 <- foreach(i=1:n_specs_1)%dopar%{
  combined_samples <- do.call("rbind",mcmc_lists_1[[i]])
  
  x_samples <- subset_x(combined_samples)%>%unlist()%>%
    array(dim=c(n_samples_1,n_obs_1,n_types))
  
  N_i <- experiment_1_specs$N[n_chains_1 * i]
  
  x_samples/N_i
}

# Relative mean samples for all locations.
mu_samples_1 <- foreach(i=1:n_specs_1)%dopar%{
  combined_samples <- do.call("rbind",mcmc_lists_1[[i]])
  
  # Coefficient samples.
  b_samples <- subset_b(combined_samples)%>%unlist()%>%
    array(dim=c(n_samples_1,experiment_1_specs$n_basis[n_chains_1*i],n_types))
  
  Z_full <- experiment_1_outputs[[n_chains_1*i]]$outputs$Z_full
  
  mu_samples <- array(NA,dim=c(n_samples_1,n_full,n_types))
  for(d in 1:n_types){
    mu_samples[,,d] <- expit(b_samples[,,d]%*%t(Z_full))
  }
  
  mu_samples
}
  
# Posterior predictive samples for all locations.
x_pred_samples_1 <- foreach(i=1:n_specs_1)%dopar%{
  combined_samples <- do.call("rbind",mcmc_lists_1[[i]])
  
  # Coefficient samples.
  b_samples <- subset_b(combined_samples)%>%unlist()%>%
    array(dim=c(n_samples_1,experiment_1_specs$n_basis[n_chains_1*i],n_types))
  
  Z_full <- experiment_1_outputs[[n_chains_1*i]]$outputs$Z_full
  
  mu_samples <- array(NA,dim=c(n_samples_1,n_full,n_types))
  for(d in 1:n_types){
    mu_samples[,,d] <- expit(b_samples[,,d]%*%t(Z_full))
  }
  
  dataset_index <- experiment_1_specs$dataset_index[n_chains_1 * i]
  train_rows <- observed_rows_1[[dataset_index]]

  ## Compute variance parameter.
  mu_mean_train <- apply(mu_samples[,train_rows,],c(1,3),mean)
  
  phi_samples <- array(NA, dim = c(n_samples_1, n_full, n_types))
  
  # Determine whether the model was run with the polynomial.
  include_polynomial <- 
    !"include_polynomial" %in% names(experiment_1_specs) ||
    isTRUE(experiment_1_specs$include_polynomial[n_chains_1*i])
  
  if(include_polynomial){
    # Polynomial coefficients.
    psi_samples <- subset_psi(combined_samples)%>%unlist()%>%
      array(dim=c(n_samples_1,4,n_types))
    
    for(d in 1:n_types){
      phi_samples[,,d] <- exp(psi_samples[,1, d] +
                                psi_samples[,2, d]*(mu_samples[,, d]-mu_mean_train[,d]) +
                                psi_samples[,3, d]*(mu_samples[,, d]-mu_mean_train[,d])^2 +
                                psi_samples[,4, d]*(mu_samples[,, d]-mu_mean_train[,d])^3)
    }
  }else{
    for(d in 1:n_types){
      phi_samples[,,d] <- combined_samples[,grepl("^phi",colnames(combined_samples))][,d]
    }
  }
  
  
  
  # Posterior predictive simulation 
  # from the Beta-Binomials.
  x_pred_samples <- array(NA, dim = c(n_samples_1, n_full, n_types))
  
  N_i <- experiment_1_specs$N[n_chains_1 * i]
  
  remainder <- matrix(N_i, nrow = n_samples_1, ncol = n_full)
  
  for(d in 1:n_types){
    p_d <- rbeta(
      n_samples_1 * n_full,
      shape1 = as.vector(mu_samples[,,d] * phi_samples[,,d]),
      shape2 = as.vector((1 - mu_samples[,,d]) * phi_samples[,,d])
    )
    
    x_d <- rbinom(
      n_samples_1 * n_full,
      size = as.vector(remainder),
      prob = p_d
    )
    
    x_d <- matrix(x_d, nrow = n_samples_1, ncol = n_full)
    
    x_pred_samples[,,d] <- x_d
    
    remainder <- remainder - x_d
  }
  
  x_pred_samples/N_i
}

# Combine direct samples from the MCMC for missing
# tree species at observed locations with posterior
# predictive samples for observed species and fully
# unobserved locations.
x_coalesced_samples_1 <- foreach(i=1:n_specs_1)%dopar%{
  
  dataset_index <- experiment_1_specs$dataset_index[n_chains_1 * i]
  
  train_rows <- sort(observed_rows_1[[dataset_index]])
  
  missing_data_i <- tree_data_with_missing_1[[dataset_index]]
  
  n_missing_i <- missing_data_i$n_missing
  
  # Start with posterior predictive samples for all rows/species.
  # This is what we want for locations with
  # n_missing == 0, locations with n_missing == 4 and
  # observed species in locations with n_missing == 1, 2, 3.
  x_coalesced <- x_pred_samples_1[[i]]
  
  # Locations with some species missing.
  partial_rows <- which(n_missing_i %in% 1:3)

  partial_train_positions <- match(partial_rows, train_rows)
  
  # Logical matrix: n_full x n_types.
  # TRUE means this species value was artificially missing.
  species_missing <- is.na(as.matrix(missing_data_i[, tree_types]))
  
  # For partially missing rows, replace only the missing species cells
  # with direct MCMC samples.
  for(k in seq_along(partial_rows)) {
    
    row_full <- partial_rows[k]
    row_train <- partial_train_positions[k]
    
    missing_species <- which(species_missing[row_full, ])
    
    if(length(missing_species) > 0) {
      x_coalesced[, row_full, missing_species] <-
        x_samples_1[[i]][, row_train, missing_species]
    }
  }
  
  x_coalesced
}

# Compute posterior predictive summaries.
pred_summaries_1 <- foreach(i = 1:n_specs_1)%dopar%{
  
  dataset_index <- experiment_1_specs$dataset_index[n_chains_1 * i]
  
  tree_cols <- c("larch", "oak", "sitka_spruce", "sycamore")
  
  # Reload the original data.
  tree_proportions <- readRDS(tree_proportions_path) %>%
    select(X_coord, Y_coord, all_of(tree_cols))
  
  missing_data_i <- tree_data_with_missing_1[[dataset_index]]
  
  # Add number of missing species column.
  tree_proportions <- tree_proportions %>%
    mutate(n_missing = missing_data_i$n_missing)
  
  # Indicators of whether each species was observed.
  observed_mat <- !is.na(as.matrix(missing_data_i[, tree_cols]))
  
  colnames(observed_mat) <- tree_cols
  
  # Original complete values as a matrix.
  actual_mat <- tree_proportions %>%
    select(all_of(tree_cols)) %>%
    as.matrix()
  
  # Use coalesced samples.
  samples_i <- x_coalesced_samples_1[[i]]
  
  # Compute CRPS values (n_full x n_types).
  crps_mat <- matrix(NA_real_, nrow = n_full, ncol = n_types)
  
  for(j in seq_len(n_full)) {
    for(d in seq_len(n_types)) {
      crps_mat[j, d] <- crps_from_samples(
        samples = samples_i[, j, d],
        y       = actual_mat[j, d]
      )
    }
  }
  
  # Posterior summaries.
  lower_mat <- apply(samples_i, c(2, 3), quantile, 0.025)
  median_mat <- apply(samples_i, c(2, 3), median)
  mean_mat <- apply(samples_i, c(2, 3), mean)
  upper_mat <- apply(samples_i, c(2, 3), quantile, 0.975)
  
  # Determine whether the model was run with the polynomial.
  include_polynomial <- 
    !"include_polynomial" %in% names(experiment_1_specs) ||
    isTRUE(experiment_1_specs$include_polynomial[n_chains_1*i])
  
  # Put all summaries together in long format.
  tree_proportions_long <- tree_proportions %>%
    pivot_longer(
      cols = all_of(tree_cols),
      names_to = "species",
      values_to = "actual"
    ) %>%
    mutate(
      observed = as.vector(t(observed_mat)),
      lower    = as.numeric(t(lower_mat)),
      median   = as.numeric(t(median_mat)),
      mean     = as.numeric(t(mean_mat)),
      upper    = as.numeric(t(upper_mat)),
      crps     = as.numeric(t(crps_mat)),
      n_basis  = experiment_1_specs$n_basis[n_chains_1 * i],
      n_obs    = experiment_1_specs$n_obs[n_chains_1 * i],
      rep      = experiment_1_specs$rep[n_chains_1 * i],
      dataset_index = dataset_index
    ) %>%
    group_by(species) %>%
    mutate(
      naive = mean(actual[observed == TRUE]),
      model = ifelse(include_polynomial,
                     "GDM+",
                     "GDM")
    ) %>%
    ungroup()
  
  tree_proportions_long
}%>%
  do.call("rbind",.)



# Quasi-Binomial GAMs -----------------------------------------------------

qb_gam_run_time <- system.time({
  qb_gam_preds <- foreach(i = 1:(n_specs_1/2)) %dopar% {
    
    # Unique experiment replicate in experiment_1_specs,
    # noting that the specs has 2 x n_chains rows
    # per replicate (due to include_polynomial.)
    fit_index <- 2 * n_chains_1 * i
    
    dataset_index <- experiment_1_specs$dataset_index[fit_index]
    
    working_data <- tree_data_with_missing_1[[dataset_index]]
    
    N_i <- experiment_1_specs$N[fit_index]
    
    n_basis_i <- experiment_1_specs$n_basis[fit_index]
    
    # Indicator matrix of whether each species was observed.
    species_observed_mat <- !is.na(as.matrix(working_data[, tree_types]))
    colnames(species_observed_mat) <- tree_types
    
    # Reload the original data.
    tree_proportions <- readRDS(tree_proportions_path) %>%
      select(X_coord, Y_coord, all_of(tree_types))
    
    actual_mat <- tree_proportions %>%
      select(all_of(tree_types)) %>%
      as.matrix()
    
    # Convert proportions to counts.
    working_data <- working_data %>%
      mutate(
        across(
          all_of(tree_types),
          ~ round(.x * N_i)
        )
      ) %>%
      mutate(
        total = pmax(
          rowSums(pick(all_of(tree_types)), na.rm = TRUE),
          N_i
        )
      )
    
    # Fit the GAMs (one per species).
    preds <- lapply(tree_types, function(sp) {
      
      form_i <- as.formula(
        paste0(
          sp,
          " / total ~ s(X_coord, Y_coord, bs = 'gp', m = c(-2, 0.5, 1), k = ",
          n_basis_i,
          ")"
        )
      )
      
      gam_i <- gam(
        form_i,
        family = quasibinomial(link = "logit"),
        weights = rep(N_i, nrow(working_data)),
        data = working_data,
        optimizer = "efs",
        method = "GCV.Cp"
      )
      
      # Predict from the GAMs.
      predict(
        gam_i,
        newdata = working_data,
        type = "response"
      )
    })
    
    pred_mat <- do.call(cbind, preds)
    colnames(pred_mat) <- tree_types
    
    # Combine all predictions and metadata
    # in a long format.
    tree_proportions_long <- tree_proportions %>%
      mutate(n_missing = working_data$n_missing) %>%
      pivot_longer(
        cols = all_of(tree_types),
        names_to = "species",
        values_to = "actual"
      ) %>%
      mutate(
        mean = as.numeric(t(pred_mat)),
        observed  = as.vector(t(species_observed_mat)),
        n_basis   = n_basis_i,
        n_obs     = experiment_1_specs$n_obs[fit_index],
        rep       = experiment_1_specs$rep[fit_index],
        crps = NA,
        lower = NA,
        upper = NA,
        dataset_index = dataset_index
      ) %>%
      group_by(species) %>%
      mutate(
        naive = mean(actual[observed == TRUE], na.rm = TRUE)
      ) %>%
      ungroup()
    
    tree_proportions_long
  }
})

qb_gam_preds%<>%do.call("rbind",.)%>%
  mutate(model="QB GAMs")

# Negative-Binomial GAMs --------------------------------------------------

nb_gam_run_time <- system.time({
  nb_gam_preds <- foreach(i = 1:(n_specs_1/2)) %dopar% {
    
    # Unique experiment replicate in experiment_1_specs,
    # noting that the specs has 2 x n_chains rows
    # per replicate (due to include_polynomial.)
    fit_index <- 2 * n_chains_1 * i
    
    dataset_index <- experiment_1_specs$dataset_index[fit_index]
    
    working_data <- tree_data_with_missing_1[[dataset_index]]
    
    N_i <- experiment_1_specs$N[fit_index]
    
    n_basis_i <- experiment_1_specs$n_basis[fit_index]
    
    # Indicator matrix of whether each species was observed.
    species_observed_mat <- !is.na(as.matrix(working_data[, tree_types]))
    colnames(species_observed_mat) <- tree_types
    
    # Reload the original data.
    tree_proportions <- readRDS(tree_proportions_path) %>%
      select(X_coord, Y_coord, all_of(tree_types))
    
    actual_mat <- tree_proportions %>%
      select(all_of(tree_types)) %>%
      as.matrix()
    
    # Convert proportions to counts.
    working_data <- working_data %>%
      mutate(
        across(
          all_of(tree_types),
          ~ round(.x * N_i )
        )
      )
    
    # Fit the GAMs (one per species).
    preds <- lapply(tree_types, function(sp) {
      
      form_i <- as.formula(
        paste0(
          sp,
          " ~ s(X_coord, Y_coord, bs = 'gp', m = c(-2, 0.5, 1), k = ",
          n_basis_i,
          ")"
        )
      )
      
      gam_i <- gam(
        form_i,
        family = "nb",
        data = working_data,
        optimizer = "efs",
        method = "GCV.Cp"
      )
      
      # Predict from the GAMs.
      predict(
        gam_i,
        newdata = working_data,
        type = "response"
      )
    })
    
    pred_mat <- do.call(cbind, preds)
    colnames(pred_mat) <- tree_types
    
    # Combine all predictions and metadata
    # in a long format.
    tree_proportions_long <- tree_proportions %>%
      mutate(n_missing = working_data$n_missing) %>%
      pivot_longer(
        cols = all_of(tree_types),
        names_to = "species",
        values_to = "actual"
      ) %>%
      mutate(
        mean = as.numeric(t(pred_mat))/N_i,
        observed  = as.vector(t(species_observed_mat)),
        n_basis   = n_basis_i,
        n_obs     = experiment_1_specs$n_obs[fit_index],
        rep       = experiment_1_specs$rep[fit_index],
        crps = NA,
        lower = NA,
        upper = NA,
        dataset_index = dataset_index
      ) %>%
      group_by(species) %>%
      mutate(
        naive = mean(actual[observed == TRUE], na.rm = TRUE)
      ) %>%
      ungroup()
    
    tree_proportions_long
  }
})

nb_gam_preds%<>%do.call("rbind",.)%>%
  mutate(model="NB GAMs")

# Beta-Binomial GAMs ------------------------------------------------------
# The BB GAMs use a lot of system memory and
# with 64Gb memory it was not possible to
# use more than two parallel workers.
registerDoParallel(cores=2)

bb_gam_run_time <- system.time({
  bb_gam_preds <- foreach(i = 1:(n_specs_1/2)) %dopar% {
    
    # Unique experiment replicate in experiment_1_specs,
    # noting that the specs has 2 x n_chains rows
    # per replicate (due to include_polynomial.)
    fit_index <- 2 * n_chains_1 * i
    
    dataset_index <- experiment_1_specs$dataset_index[fit_index]
    
    working_data <- tree_data_with_missing_1[[dataset_index]]
    
    N_i <- experiment_1_specs$N[fit_index]
    
    n_basis_i <- experiment_1_specs$n_basis[fit_index]
    
    # Indicator matrix of whether each species was observed.
    species_observed_mat <- !is.na(as.matrix(working_data[, tree_types]))
    colnames(species_observed_mat) <- tree_types
    
    # Reload the original data.
    tree_proportions <- readRDS(tree_proportions_path) %>%
      select(X_coord, Y_coord, all_of(tree_types))
    
    actual_mat <- tree_proportions %>%
      select(all_of(tree_types)) %>%
      as.matrix()
    
    # Convert proportions to counts.
    # Convert proportions to counts.
    working_data <- working_data %>%
      mutate(
        across(
          all_of(tree_types),
          ~ round(.x * N_i)
        )
      ) %>%
      mutate(
        total = pmax(
          rowSums(pick(all_of(tree_types)), na.rm = TRUE),
          N_i
        )
      )
    
    # Fit the GAMs (one per species).
    preds <- lapply(tree_types, function(sp) {
      
      form_i <- as.formula(
        paste0(
          "cbind(",
          sp,
          ", total - ",
          sp,
          ") ~ s(X_coord, Y_coord, bs = 'gp', m = c(-2, 0.5, 1), k = ",
          n_basis_i,
          ")"
        )
      )
      
      # Fit the GAM using glmmTMB.
      gam_i <-  glmmTMB(
        form_i,
        family = betabinomial(link = "logit"),
        data = working_data,
        REML = TRUE
      )
      
      # Predict from the GAMs.
      bb_preds <- predict(
        gam_i,
        newdata = working_data,
        type = "link",
        se.fit = FALSE,
        cov.fit = FALSE
      )%>%
        plogis()
      
      # Delete the current GAM to
      # reduce memory usage.
      rm(gam_i)
      gc()
      
      bb_preds
    })
    
    pred_mat <- do.call(cbind, preds)
    colnames(pred_mat) <- tree_types
    
    # Combine all prediction summaries
    # together in a long format.
    tree_proportions_long <- tree_proportions %>%
      mutate(n_missing = working_data$n_missing) %>%
      pivot_longer(
        cols = all_of(tree_types),
        names_to = "species",
        values_to = "actual"
      ) %>%
      mutate(
        mean = as.numeric(t(pred_mat)),
        observed  = as.vector(t(species_observed_mat)),
        n_basis   = n_basis_i,
        n_obs     = experiment_1_specs$n_obs[fit_index],
        rep       = experiment_1_specs$rep[fit_index],
        crps = NA,
        lower = NA,
        upper = NA,
        dataset_index = dataset_index
      ) %>%
      group_by(species) %>%
      mutate(
        naive = mean(actual[observed == TRUE], na.rm = TRUE)
      ) %>%
      ungroup()
    
    tree_proportions_long
  }
})

bb_gam_preds%<>%do.call("rbind",.)%>%
  mutate(model="BB GAMs")


# Compare prediction performance ------------------------------------------

# Compute performance summaries.
compare_preds <- pred_summaries_1%>%
  full_join(qb_gam_preds)%>%
  full_join(nb_gam_preds)%>%
  full_join(bb_gam_preds)%>%
  mutate(model=factor(model,
                      levels=c("GDM+","GDM","BB GAMs",
                               "QB GAMs","NB GAMs")))

pred_perf_1 <- compare_preds%>%
  filter(observed==FALSE)%>%
  group_by(n_missing,model)%>%
  summarise(mae=mean(abs(actual-mean)),
            rmse=sqrt(mean((actual-mean)^2)),
            mean_crps=mean(crps),
            coverage=mean(actual>=lower&actual<=upper),
            mean_width=mean(upper-lower),
            xi=1-mean((actual-mean)^2)/mean((actual-naive)^2),
            n=n())

pred_perf_1 

# Compute performance summaries by species.
pred_perf_by_species_1 <- compare_preds%>%
  filter(observed==FALSE)%>%
  group_by(n_missing,model,species)%>%
  summarise(mae=mean(abs(actual-mean)),
            rmse=sqrt(mean((actual-mean)^2)),
            mean_crps=mean(crps),
            coverage=mean(actual>=lower&actual<=upper),
            mean_width=mean(upper-lower),
            xi=1-mean((actual-mean)^2)/mean((actual-naive)^2),
            n=n())

# Create full table for the paper.
table_experiment_1 <- full_join(pred_perf_1%>%mutate(species="Overall"),
                                pred_perf_by_species_1)%>%
  mutate(species=case_match(species,"larch"~"Larch",
                         "oak"~"Oak","sitka_spruce"~"Sitka spruce",
                         "sycamore"~"Sycamore",.default = species),
  species=factor(species,levels=c("Larch","Oak","Sitka spruce","Sycamore","Overall")))%>%
  arrange(species,n_missing)%>%
  select(species,n_missing,n,model,mae,rmse,xi,mean_crps,mean_width,coverage)

write.xlsx(table_experiment_1,file="Outputs/experiment_1_results_full.xlsx",asTable = TRUE)

table_compare_models <- table_experiment_1%>%
  filter(model%in%c("GDM+","BB GAMs"))%>%
  select(-mean_crps,-coverage,-mean_width)%>%
  pivot_wider(names_from = model,values_from = c(mae,rmse,xi))

write.xlsx(table_compare_models,file="Outputs/table_compare_models.xlsx",asTable = TRUE)

# Scatter plot for an example replicate.
scatter_compare <- compare_preds %>%
  filter(model %in% c("GDM+","BB GAMs"),
         n_missing<=3,observed==FALSE,rep==1) %>%
  mutate(n_missing=case_match(n_missing,
                              1~"1 species missing",
                              2~"2 species missing",
                              3~"3 species missing"),
         species = case_when(species == 'larch' ~ 'Larch',
                          species == 'oak' ~ 'Oak',
                          species == 'sitka_spruce' ~ 'Sitka spruce',
                          species == 'sycamore' ~ 'Sycamore',
                          TRUE ~ species)) %>%
  
  ggplot(., aes(x = mean, y = actual,
                shape=species)) +
  scale_x_continuous(limits=c(0,1))+
  scale_y_continuous(limits=c(0,1))+
  geom_abline() +
  geom_point() +
  labs(x="Predicted proportion",y="Observed proportion",
       shape=NULL,colour=NULL)+
  facet_grid(model ~ n_missing)+
  #coord_fixed()+
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(vjust=0.8))+
  scale_shape_manual(values=c(c(0,1,2,5)))+
  scale_colour_viridis_d(begin = 0.1,end=0.9)
ggsave(scatter_compare,filename="Plots/scatter_compare.pdf",width=8,height=5)

# Line plot of MAE and RMSE by model.
mae_rmse_plot <- ggplot(pred_perf_1%>%
                     #filter(model!="NB GAMs")%>%
                     select(model,n_missing,rmse,mae)%>%
                     rename(`Root mean square error`=rmse,
                            `Mean absolute error`=mae)%>%
                     pivot_longer(cols=c("Mean absolute error","Root mean square error"),
                                  names_to = "metric")%>%
                     mutate(model=factor(model,levels=c("NB GAMs","QB GAMs","BB GAMs","GDM","GDM+"))),
                   aes(x=n_missing,y=value,colour = model))+
  facet_wrap(~metric,scales="free")+
  geom_line()+
  geom_point()+
  theme(panel.grid = element_blank())+
  labs(x="Number of species missing",y=NULL,
       colour="Method")+
  scale_colour_viridis_d(begin = 0.1,end=0.9,direction=-1)
ggsave(mae_rmse_plot,filename="Plots/mae_rmse_plot.pdf",width=8,height=3)

# Check zero proportions --------------------------------------------------

# Context lookup.
missing_context_lookup <- tibble(
  n_missing = c(1, 2, 3, 4),
  context = c(
    "1 species missing",
    "2 species missing",
    "3 species missing",
    "All species missing"
  )
)

n_missing_labels_eval <- missing_context_lookup$context

# Compute posterior predictive/imputed proportion of zeros
# for focal species values that were actually made missing.
zero_props_1 <- foreach(i = seq_len(n_specs_1), .combine = bind_rows) %dopar% {
  
  dataset_index <- experiment_1_specs$dataset_index[n_chains_1 * i]
  
  missing_data_i <- tree_data_with_missing_1[[dataset_index]]
  n_missing_i <- missing_data_i$n_missing
  
  # Matrix identifying whether each species at each location 
  # has a missing value or not.
  species_missing_mat <- is.na(as.matrix(missing_data_i[, tree_types]))
  colnames(species_missing_mat) <- tree_types
  
  samples_i <- x_coalesced_samples_1[[i]]
  
  # Determine whether the model was run with the polynomial.
  include_polynomial <- 
    !"include_polynomial" %in% names(experiment_1_specs) ||
    isTRUE(experiment_1_specs$include_polynomial[n_chains_1 * i])
  
  
  zero_prop_array <- array(dim=c(n_samples_1,n_types,4))
  
  # For each MCMC sample.
  for(ii in 1:n_samples_1){
    # For each species.
    for(j in 1:n_types){
      # For each number of missing species.
      for(k in 1:4){
        # Identify rows/locations where this species is missing
        # and where the number of missing species is k.
        eval_rows <- n_missing_i==k & species_missing_mat[,j]==TRUE
        
        # Compute zero proportion in this subset.
        zero_prop_array[ii,j,k] <- mean(samples_i[ii,eval_rows,j]==0)
      }
    }
  }
  
  # Make into a data frame.
  zero_prop_df <- zero_prop_array%>%melt()%>%
    mutate(sample=Var1,species=tree_types[Var2],n_missing=Var3,zero_prop=value)%>%
    left_join(missing_context_lookup, by = "n_missing") %>%
    mutate(species_index = match(species, tree_types),
           n_basis = experiment_1_specs$n_basis[n_chains_1 * i],
           n_obs = experiment_1_specs$n_obs[n_chains_1 * i],
           rep = experiment_1_specs$rep[n_chains_1 * i],
           dataset_index = dataset_index,
           model = ifelse(include_polynomial, "GDM+", "GDM")
    ) %>%
    select(
      sample,
      species,
      n_missing,
      context,
      zero_prop,
      n_basis,
      n_obs,
      rep,
      dataset_index,
      model
    )
  
  zero_prop_df
}

# Reload the original complete data.
truth_full <- readRDS(tree_proportions_path) %>%
  select(all_of(tree_types)) %>%
  as.matrix()

# Compute true proportion of zeros for the same focal species cells
# that were artificially made missing.
zero_props_truth_1 <- foreach(i = seq_len(n_specs_1), .combine = bind_rows) %dopar% {
  
  dataset_index <- experiment_1_specs$dataset_index[n_chains_1 * i]
  
  missing_data_i <- tree_data_with_missing_1[[dataset_index]]
  n_missing_i <- missing_data_i$n_missing
  
  # Matrix identifying whether each species at each location 
  # has a missing value or not.
  species_missing_mat <- is.na(as.matrix(missing_data_i[, tree_types]))
  colnames(species_missing_mat) <- tree_types
  
  # Determine whether the model was run with the polynomial.
  include_polynomial <- 
    !"include_polynomial" %in% names(experiment_1_specs) ||
    isTRUE(experiment_1_specs$include_polynomial[n_chains_1 * i])
  
  zero_prop_array <- array(dim=c(n_types,4))
  
  for(j in 1:n_types){
    for(k in 1:4){
      # Identify rows/locations where this species is missing
      # and where the number of missing species is k.
      eval_rows <- n_missing_i==k & species_missing_mat[,j]==TRUE
      
      # Compute zero proportion in this subset.
      zero_prop_array[j,k] <- mean(truth_full[eval_rows,tree_types[j]]==0)
    }
  }
  
  # Make into a data frame.
  zero_prop_df <- zero_prop_array%>%melt()%>%
    mutate(species=tree_types[Var1],n_missing=Var2,zero_prop=value)%>%
    left_join(missing_context_lookup, by = "n_missing") %>%
    mutate(species_index = match(species, tree_types),
           n_basis = experiment_1_specs$n_basis[n_chains_1 * i],
           n_obs = experiment_1_specs$n_obs[n_chains_1 * i],
           rep = experiment_1_specs$rep[n_chains_1 * i],
           dataset_index = dataset_index,
           model = ifelse(include_polynomial, "GDM+", "GDM")
    ) %>%
    select(
      species,
      n_missing,
      context,
      zero_prop,
      n_basis,
      n_obs,
      rep,
      dataset_index,
      model
    )
  
  zero_prop_df
}

# Summarise the posterior predictive distributions with 
# means and 95% intervals.
zero_props_summary_1 <- zero_props_1 %>%
  group_by(model, species, n_missing, context, rep, dataset_index) %>%
  summarise(
    lower = quantile(zero_prop, 0.025, na.rm = TRUE),
    mean = mean(zero_prop, na.rm = TRUE),
    upper = quantile(zero_prop, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    zero_props_truth_1,
    by = c(
      "species",
      "n_missing",
      "context",
      "rep",
      "dataset_index",
      "model"
    ),
    suffix = c("", "_truth")
  ) %>%
  mutate(
    context = factor(context, levels = n_missing_labels_eval),
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
zero_prop_plot <- ggplot(zero_props_summary_1%>%
                           filter(model=="GDM+"), aes(x = mean, y = zero_prop, 
                                 shape = species)) +
  facet_wrap(~ context) +
  geom_abline() +
  geom_errorbar(aes(xmin = lower, xmax = upper)) +
  geom_point() +
  theme(panel.grid = element_blank())+
  labs(
    x = "Predicted zero proportion",
    y = "Observed zero proportion",
    colour=NULL,shape=NULL
  ) +
  scale_shape_manual(values=c(0,1,2,5))+
  scale_x_continuous(limits=c(0.3,0.9))+
  scale_y_continuous(limits=c(0.3,0.9))
ggsave(zero_prop_plot,filename="Plots/zero_prop.pdf",width=8,height=5)

# Overall zero proportion 95% P.I. coverage.
zero_props_summary_1%>%
  group_by(model)%>%
  summarise(mean(zero_prop>=lower&zero_prop<=upper))

# Zero proportion MAE by number
# of missing species.
zero_props_summary_1%>%
  group_by(context,model)%>%
  summarise(mean(abs(zero_prop-mean)))

# Zero proportion MAE by number of missing
# species and by species.
zero_props_summary_1%>%
  group_by(context,species,model)%>%
  summarise(mean(abs(zero_prop-mean)))

# Zero proportion MAE by species.
zero_props_summary_1%>%
  group_by(species,model)%>%
  summarise(mean(abs(zero_prop-mean)))

# Overall MAE, RMSE, mean CRPS for
# the GDM+ and GDM models.
compare_GDM_overall <- compare_preds%>%
  filter(model%in%c("GDM","GDM+"))%>%
  group_by(model)%>%
  summarise(mae=mean(abs(actual-mean)),
            rmse=sqrt(mean((actual-mean)^2)),
            mean_crps=mean(crps))

# Compare RMSE and MAE.
compare_GDM_overall[1,2]/compare_GDM_overall[2,2]
1-compare_GDM_overall[1,c(-1,-2)]/compare_GDM_overall[2,c(-1,-2)]

# Compute zero proportion bias, MAE, and RMSE.
compare_zero_overall <- zero_props_summary_1%>%
  group_by(model)%>%
  summarise(zero_bias=mean(zero_prop-mean),
            zero_mae=mean(abs(zero_prop-mean)),
            zero_rmse=sqrt(mean((zero_prop-mean)^2)))

# Ratio between models.
1-compare_zero_overall[2,c(2,3,4)]/compare_zero_overall[1,c(2,3,4)]

# Zero proportion bias, MAE, and RMSE
# by number of missing species.
zp_perf_1 <- zero_props_summary_1%>%
  group_by(model,n_missing)%>%
  summarise(zero_bias=mean(zero_prop-mean),
            zero_mae=mean(abs(zero_prop-mean)),
             zero_rmse=sqrt(mean((zero_prop-mean)^2)))

# Table comparing GDM performance.
table_compare_GDMs <- table_experiment_1%>%
  filter(species=="Overall")%>%
  filter(model%in%c("GDM+","GDM"))%>%
  select(-mae,-rmse,-xi,-mean_width,-n,-species)%>%
  left_join(select(zp_perf_1,-zero_rmse))%>%
  pivot_wider(names_from = model,values_from = c(mean_crps,coverage,zero_bias,zero_mae))

write.xlsx(table_compare_GDMs,file="Outputs/table_compare_GDMs.xlsx",asTable = TRUE)

# Example heatmap ---------------------------------------------------------

# Convert Beta-Binomial mean parameter values
# from repeat 1 to absolute means for
# each species.
mean_samples_1 <- mu_samples_1[[1]]%>%apply(c(1,2),function(u){
  v <- numeric(length(u))
  v[1] <- u[1]
  for(j in 2:length(u)) v[j] <- u[j]*(1-sum(v[1:(j-1)]))
  v
})%>%aperm(c(2,3,1))

mean_samples_df_1 <- mean_samples_1%>%apply(c(2,3),mean)
dimnames(mean_samples_df_1) <- list(location=1:n_full,species=tree_types)

mean_samples_df_1%<>%reshape2::melt(value.name = "mean")%>%
  mutate(model="GDM+")%>%
  cbind(select(tree_proportions,X_coord,Y_coord))%>%
  select(-location)

# Put original proportions into a long format.
tree_proportions_long <- tree_proportions %>%
  pivot_longer(
    cols = all_of(tree_types),
    names_to = "species",
    values_to = "mean"
  )%>%
  mutate(model="Original data")

# Plot heatmap of original data alongside 
# predictions from the GDM+ and BB GAMs.
mean_preds_1 <- full_join(mean_samples_df_1,
                          select(bb_gam_preds,X_coord,Y_coord,species,mean,model))%>%
  full_join(tree_proportions_long)

png("Plots/heatmap_predictions_means-all.png", height = 8,width=8,units="in",res=600)
plot(mean_preds_1 %>%
       mutate(species = case_when(species == 'larch' ~ 'Larch',
                               species == 'oak' ~ 'Oak',
                               species == 'sitka_spruce' ~ 'Sitka spruce',
                               species == 'sycamore' ~ 'Sycamore',
                               TRUE ~ species)) %>%
       mutate(model = factor(model, levels = c('Original data', 'GDM+', 'BB GAMs'))) %>%
       
       ggplot(., aes(x = X_coord, y = Y_coord, fill = mean)) +
       geom_tile() +
       scale_fill_viridis_c(name="Proportion",breaks=seq(0,1,by=0.2)) +scale_colour_viridis_c(name="Proportion",breaks=seq(0,1,by=0.2)) +
       theme(axis.text.x=element_text(angle=60, hjust=1)) +
       coord_fixed()+
       facet_grid(species ~ model) +
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
              legend.title = element_text(vjust=0.8)) 
)
dev.off()


# Polynomial effect plot --------------------------------------------------

# Posterior predictive samples for all locations.
concentration_parameter_samples_1 <- foreach(i=1:n_specs_1)%dopar%{
  combined_samples <- do.call("rbind",mcmc_lists_1[[i]])
  
  mu_seq <- seq(0.001,0.999,by=0.001)
  n_mu <- length(mu_seq)
  
  # Coefficient samples.
  b_samples <- subset_b(combined_samples)%>%unlist()%>%
    array(dim=c(n_samples_1,experiment_1_specs$n_basis[n_chains_1*i],n_types))
  
  Z_full <- experiment_1_outputs[[n_chains_1*i]]$outputs$Z_full
  
  mu_samples <- array(NA,dim=c(n_samples_1,n_full,n_types))
  for(d in 1:n_types){
    mu_samples[,,d] <- expit(b_samples[,,d]%*%t(Z_full))
  }
  
  dataset_index <- experiment_1_specs$dataset_index[n_chains_1 * i]
  train_rows <- observed_rows_1[[dataset_index]]
  
  ## Compute variance parameter.
  mu_mean_train <- apply(mu_samples[,train_rows,],c(1,3),mean)
  
  phi_samples <- array(NA, dim = c(n_samples_1, n_mu, n_types))
  
  # Determine whether the model was run with the polynomial.
  include_polynomial <- 
    !"include_polynomial" %in% names(experiment_1_specs) ||
    isTRUE(experiment_1_specs$include_polynomial[n_chains_1*i])
  
  if(include_polynomial){
    # Polynomial coefficients.
    psi_samples <- subset_psi(combined_samples)%>%unlist()%>%
      array(dim=c(n_samples_1,4,n_types))
    
    for (d in seq_len(n_types)) {
      
      mu_centred <- matrix(
        rep(mu_seq, each = n_samples_1),
        nrow = n_samples_1,
        ncol = n_mu
      ) - mu_mean_train[, d]
      
      phi_samples[, , d] <- exp(
        psi_samples[, 1, d] +
          psi_samples[, 2, d] * mu_centred +
          psi_samples[, 3, d] * mu_centred^2 +
          psi_samples[, 4, d] * mu_centred^3
      )
    }
  }else{
    for(d in 1:n_types){
      phi_samples[,,d] <- combined_samples[,grepl("^phi",colnames(combined_samples))][,d]
    }
  }

  phi_quantiles_df <- phi_samples%>%apply(c(2,3),quantile,c(0.025,0.5,0.975))%>%
    melt()%>%pivot_wider(names_from="Var1",values_from="value")%>%
    rename(lower=`2.5%`,median=`50%`,upper=`97.5%`,
           mu=Var2,species=Var3)%>%
    mutate(model = ifelse(include_polynomial,
                          "GDM+",
                          "GDM"),
           rep=experiment_1_specs$rep[n_chains_1*i],
           species=tree_types[species],
           mu=mu_seq[mu])
  
  phi_quantiles_df
}%>%
  do.call("rbind",.)

phi_plot <- ggplot(concentration_parameter_samples_1%>%
         filter(rep==1)%>%
         mutate(species = case_when(
           species == "larch" ~ "Larch",
           species == "oak" ~ "Oak",
           species == "sitka_spruce" ~ "Sitka spruce",
           species == "sycamore" ~ "Sycamore",
           TRUE ~ species
         )),
       aes(x=mu,y=median,ymin=lower,ymax=upper))+
  facet_wrap(~species,nrow=1)+
  geom_ribbon(aes(fill=model),alpha=0.5)+
  geom_line(aes(colour=model))+
  scale_x_continuous(limits=c(-0.05,1.05))+
  scale_y_log10()+
  theme(panel.grid = element_blank())+
  labs(
    x = expression("Beta-Binomial mean parameter "*mu["s,d"]),
    y = expression("Concentration parameter "*phi["s,d"]),
    colour=NULL,shape=NULL,fill=NULL
  ) +
  scale_colour_viridis_d(begin = 0.1,end=0.9)+
  scale_fill_viridis_d(begin = 0.1,end=0.9)
ggsave(phi_plot,file="Plots/phi_plot.pdf",width=8,height=2.5)

# Save outputs ------------------------------------------------------------

# Save pretty much everything.
save(experiment_1_specs,observed_rows_1,
     n_specs_1,n_chains_1,n_samples_1,n_obs_1,
     total_run_time_1,experiment_1_outputs,
     tree_data_with_missing_1,
     x_pred_samples_1,pred_summaries_1,
     x_coalesced_samples_1,mu_samples_1,
     zero_props_1,zero_props_truth_1,
     zero_props_summary_1,
     qb_gam_preds,qb_gam_run_time,
     nb_gam_preds,nb_gam_run_time,
     bb_gam_preds,bb_gam_run_time,
     compare_preds,
     pred_perf_1,file="Saved runs/experiment_1_10_06_26.RData")

# Save only smaller objects required to inspect
# and plot the results.
save(experiment_1_specs,observed_rows_1,
     n_specs_1,n_chains_1,n_samples_1,n_obs_1,
     total_run_time_1,
     tree_data_with_missing_1,
     zero_props_summary_1,
     mean_preds_1,
     qb_gam_preds,qb_gam_run_time,
     nb_gam_preds,nb_gam_run_time,
     bb_gam_preds,bb_gam_run_time,
     compare_preds,
     pred_perf_1,file="Saved runs/experiment_1_lite_10_06_26.RData")
