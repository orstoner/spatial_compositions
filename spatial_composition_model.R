# Beta-Binomial NIMBLE functions ------------------------------------------

# Probability mass function.
dbetabinomial = nimbleFunction(
  
  run = function(x = double(0),
                 nu = double(0),
                 phi = double(0),
                 size = double(0),
                 log = integer(0)){  
    
    returnType(double(0))
    
    if(x >= 0 & x <= size){  
      
      return(lgamma(size + 1) + lgamma(x + nu * phi) + lgamma(size - x + (1 - nu) * phi) + lgamma(phi) -                
               lgamma(size + phi) - lgamma(nu * phi) - lgamma((1 - nu) * phi) - lgamma(size - x + 1) - lgamma(x + 1))     
      
    }else{       
      return(-Inf)     
    }   
  },
  buildDerivs = TRUE)     

# Simulation function.
rbetabinomial = nimbleFunction(
  
  run = function(n = integer(0),
                 nu = double(0),
                 phi = double(0),
                 size = double(0)){   
    
    returnType(double(0)) 
    
    # Upper limit on phi for computational stability.
    phi <- min(phi, 1e+04) 
    pi = rbeta(1, nu * phi, (1 - nu) * phi)     
    
    return(rbinom(1, size, pi))   
    
  })      


# Register as a distribution for NIMBLE.
registerDistributions(
  list(dbetabinomial = list(BUGSdist='dbetabinomial(nu, phi, size)', discrete=TRUE))
)

# Spatial composition model function --------------------------------------

spatial_composition_model <- function(X,
                     initial_seed, chain_seed,
                     counts_full,
                     locations_full,
                     totals_full, n_types, n_basis,
                     niter, nburn, nthin,
                     define_smooth_using_train=TRUE) {
  
  library(nimble)

  # Identify which locations have at
  # at least one observed species.
  obs_index <- which(apply(counts_full,1,function(u)any(!is.na(u))))
  
  # Create training data frames (i.e. exclude
  # locations with no data).
  n_train <- length(obs_index)
  counts_train <- counts_full[obs_index,]
  locations_train <- locations_full[obs_index,]
  totals_train <- totals_full[obs_index]

  # Set up spline basis matrix.
  sm <- mgcv::smoothCon(
    mgcv::s(X_coord,Y_coord, bs = "gp", m=c(-2,0.5,1), k = n_basis),
    data = {if(define_smooth_using_train){locations_train}else{locations_full}},
    absorb.cons = TRUE,
    diagonal.penalty = TRUE
  )[[1]]
  
  # Basis matrix over the full spatial grid.
  Z_full  <- cbind(1,mgcv::PredictMat(sm, locations_full))
  
  # Basis matrix over only the training locations.
  Z_train <- Z_full[obs_index,]
  
  # Manually center basis functions if not defined using
  # only training locations.
  if(!define_smooth_using_train){
    train_col_means <- colMeans(Z_train)
    
    Z_full <- sweep(Z_full, 2, train_col_means, "-")
    Z_full[,1] <- 1
    
    Z_train <- Z_full[obs_index, , drop = FALSE]
    
  }

  # NIMBLE code.
  spatial_model_code <- nimbleCode({
    
    # GDM likelihood.
    for (i in 1:n) {
      # First tree species.
      x[i, 1] ~ dbetabinomial(mu[i, 1], phi[i, 1], N[i])
      
      # Subsequent species.
      for (k in 2:n_types) {
        x[i, k] ~ dbetabinomial(mu[i, k], phi[i, k], N[i] - sum(x[i, 1:(k-1)]))
      }
    }

    # Loop over tree species.
    for (k in 1:n_types) {
      
      # Beta-Binomial means.
      mu[1:n, k] <- expit(eta[1:n, k])
      
      # Linear predictors.
      eta[1:n, k] <- Z[1:n, 1:n_basis] %*% b[1:n_basis, k]
      
      # Spatial splines.
      b[1,k] ~ dnorm(0, sd = 10)
      for(i in 2:n_basis){
        b[i,k] ~ dnorm(0,sd=sigma[k])
      }
      
      # Spline variance parameter priors.
      sigma[k] ~ T(dnorm(0,sd=10),0,)
      
      # Beta-Binomial inverse-dispersion parameters.
      mu_mean[k] <- mean(mu[1:n, k])
      
      # Cubic polynomials.
      for (i in 1:n) {
        log(phi[i, k]) <- psi[1, k] +
          psi[2, k]*(mu[i, k]-mu_mean[k]) +
          psi[3, k]*(mu[i, k]-mu_mean[k])^2 +
          psi[4, k]*(mu[i, k]-mu_mean[k])^3
      }
      
      # Beta-Binomial inverse-dispersion intercept prior.
      psi[1, k] ~ dnorm(2,sd=2)
      
      # Inverse-dispersion coefficient priors.
      for(j in 2:4){
        psi[j, k] ~ dnorm(0,sd=1)
      }
      
    }
    
  })
  
  # NIMBLE data.
  data <- list(x = counts_full[obs_index,])

  # NIMBLE constants.
  constants <- list(Z = Z_train,
                    n = n_train,
                    N = totals_train,
                    n_basis = n_basis,
                    n_types = n_types)
  
  
  # Function to generate initial values.
  init_fun <- function(){
    x_inits <- matrix(NA, nrow = n_train, ncol = n_types) 
    for(i in 1:n_train){
      # Initial values for missing tree species.
      if(is.na(sum(counts_train[i,]))){
        n_missing <- sum(is.na(counts_train[i,]))
        x_inits[i,is.na(counts_train[i,])] <- rmulti(1,totals_train[i]-sum(counts_train[i,which(!is.na(counts_train[i,]))]),
                                                     prob=rep(1/(n_missing+1),n_missing+1))[1:n_missing]
      }
    }
    
    inits <- list(b = matrix(rnorm(n_basis*n_types, 0, 0.05), nrow = n_basis, ncol = n_types),
                  sigma = runif(n_types,0,10),
                  psi = matrix(rnorm(4*n_types,rep(c(2,0,0,0),n_types),0.25), nrow = 4, ncol = n_types),
                  x = x_inits)
    return(inits)
  }
  
  set.seed(initial_seed[X])
  initial_values <- init_fun()
  
  # Build the NIMBLE model.
  model <- nimbleModel(code = spatial_model_code,
                       data = data,
                       constants = constants,
                       inits = init_fun())
  
  # Compile the model.
  compiled_model <- compileNimble(model)
  
  # Configure the MCMC.
  conf_model <- configureMCMC(compiled_model, 
                              monitors = c("x",
                                           "b",
                                           "sigma",
                                           "eta",
                                           "mu",
                                           "phi",
                                           "psi"
                              ),
                              print = TRUE, 
                              useConjugacy = FALSE)
  
  
  # Replace samplers for sigma and psi.
  conf_model$removeSamplers(c("sigma","psi"))

  for(j in 1:n_types){
    conf_model$addSampler(target = paste0("psi[1:",4,",",j,"]"), type = "AF_slice")
  }
  for (node in 1:length(model$expandNodeNames("sigma"))) {
    conf_model$addSampler(target = model$expandNodeNames("sigma")[node], type = "RW",control=list(log=TRUE))
  }
  
  # Build the MCMC.
  spatial_model <- buildMCMC(conf_model)
  
  # Compile the MCMC.
  compiled_spatial_model <- compileNimble(spatial_model,
                                          project = model,
                                          resetFunctions = TRUE)
  
  # Run the MCMC.
  samples <- runMCMC(compiled_spatial_model, 
                     niter =   niter, 
                     nburnin = nburn,
                     thin =   nthin, 
                     samplesAsCodaMCMC = TRUE,
                     setSeed = chain_seed[X])
  
  # Return outputs.
  return(list(initial_values = initial_values, samples = samples, Z_full=Z_full, Z_train=Z_train))
  
}


# Spatial composition model function --------------------------------------
# (no polynomial)
spatial_composition_model_no_polynomial <- function(X,
                                                    initial_seed, chain_seed,
                                                    counts_full,
                                                    locations_full,
                                                    totals_full, n_types, n_basis,
                                                    niter, nburn, nthin,
                                                    define_smooth_using_train=TRUE) {
  
  library(nimble)
  
  # Identify which locations have at
  # at least one observed species.
  obs_index <- which(apply(counts_full,1,function(u)any(!is.na(u))))
  
  # Create training data frames (i.e. exclude
  # locations with no data).
  n_train <- length(obs_index)
  counts_train <- counts_full[obs_index,]
  locations_train <- locations_full[obs_index,]
  totals_train <- totals_full[obs_index]
  
  # Set up spline basis matrix.
  sm <- mgcv::smoothCon(
    mgcv::s(X_coord,Y_coord, bs = "gp", m=c(-2,0.5,1), k = n_basis),
    data = {if(define_smooth_using_train){locations_train}else{locations_full}},
    absorb.cons = TRUE,
    diagonal.penalty = TRUE
  )[[1]]
  
  # Basis matrix over the full spatial grid.
  Z_full  <- cbind(1,mgcv::PredictMat(sm, locations_full))
  
  # Basis matrix over only the training locations.
  Z_train <- Z_full[obs_index,]
  
  # Manually center basis functions if not defined using
  # only training locations.
  if(!define_smooth_using_train){
    train_col_means <- colMeans(Z_train)
    
    Z_full <- sweep(Z_full, 2, train_col_means, "-")
    Z_full[,1] <- 1
    
    Z_train <- Z_full[obs_index, , drop = FALSE]
    
  }
  
  # NIMBLE code.
  spatial_model_code <- nimbleCode({
    
    # GDM likelihood.
    for (i in 1:n) {
      # First tree species.
      x[i, 1] ~ dbetabinomial(mu[i, 1], phi[1], N[i])
      
      # Subsequent species.
      for (k in 2:n_types) {
        x[i, k] ~ dbetabinomial(mu[i, k], phi[k], N[i] - sum(x[i, 1:(k-1)]))
      }
    }
    
    # Loop over tree species.
    for (k in 1:n_types) {
      
      # Beta-Binomial means.
      mu[1:n, k] <- expit(eta[1:n, k])
      
      # Linear predictors.
      eta[1:n, k] <- Z[1:n, 1:n_basis] %*% b[1:n_basis, k]
      
      # Spatial splines.
      b[1,k] ~ dnorm(0, sd = 10)
      for(i in 2:n_basis){
        b[i,k] ~ dnorm(0,sd=sigma[k])
      }
      
      # Spline variance parameter priors.
      sigma[k] ~ T(dnorm(0,sd=10),0,)
      
      # Beta-Binomial inverse-dispersion intercept prior.
      psi[k] ~ dnorm(2,sd=2)
      log(phi[k]) <- psi[k]
      
    }
    
  })
  
  # NIMBLE data.
  data <- list(x = counts_full[obs_index,])
  
  # NIMBLE constants.
  constants <- list(Z = Z_train,
                    n = n_train,
                    N = totals_train,
                    n_basis = n_basis,
                    n_types = n_types)
  
  
  # Function to generate initial values.
  init_fun <- function(){
    x_inits <- matrix(NA, nrow = n_train, ncol = n_types) 
    for(i in 1:n_train){
      # Initial values for missing tree species.
      if(is.na(sum(counts_train[i,]))){
        n_missing <- sum(is.na(counts_train[i,]))
        x_inits[i,is.na(counts_train[i,])] <- rmulti(1,totals_train[i]-sum(counts_train[i,which(!is.na(counts_train[i,]))]),
                                                     prob=rep(1/(n_missing+1),n_missing+1))[1:n_missing]
      }
    }
    
    inits <- list(b = matrix(rnorm(n_basis*n_types, 0, 0.05), nrow = n_basis, ncol = n_types),
                  sigma = runif(n_types,0,10),
                  psi = rnorm(n_types,2,0.25),
                  x = x_inits)
    return(inits)
  }
  
  set.seed(initial_seed[X])
  initial_values <- init_fun()
  
  # Build the NIMBLE model.
  model <- nimbleModel(code = spatial_model_code,
                       data = data,
                       constants = constants,
                       inits = init_fun())
  
  
  # Compile the model.
  compiled_model <- compileNimble(model)
  
  # Configure the MCMC.
  conf_model <- configureMCMC(compiled_model, 
                              monitors = c("x",
                                           "b",
                                           "sigma",
                                           "eta",
                                           "mu",
                                           "phi",
                                           "psi"
                              ),
                              print = TRUE, 
                              useConjugacy = FALSE)
  
  
  # Replace samplers for sigma and psi.
  conf_model$removeSamplers(c("sigma"))
  
  for (node in 1:length(model$expandNodeNames("sigma"))) {
    conf_model$addSampler(target = model$expandNodeNames("sigma")[node], type = "RW",control=list(log=TRUE))
  }
  
  # Build the MCMC.
  spatial_model <- buildMCMC(conf_model)
  
  # Compile the MCMC.
  compiled_spatial_model <- compileNimble(spatial_model,
                                          project = model,
                                          resetFunctions = TRUE)
  
  # Run the MCMC.
  samples <- runMCMC(compiled_spatial_model, 
                     niter =   niter, 
                     nburnin = nburn,
                     thin =   nthin, 
                     samplesAsCodaMCMC = TRUE,
                     setSeed = chain_seed[X])
  
  # Return outputs.
  return(list(initial_values = initial_values, samples = samples, Z_full=Z_full, Z_train=Z_train))
  
}



# Function to extract parameter samples -----------------------------------

subset_b <- function(chain) {
  cols <- grep("^b", colnames(chain), value = TRUE)
  return(chain[, cols])
}

subset_lambda <- function(chain) {
  cols <- grep("^lambda", colnames(chain), value = TRUE)
  return(chain[, cols])
}

subset_phi <- function(chain) {
  cols <- grep("^phi", colnames(chain), value = TRUE)
  return(chain[, cols])
}

subset_sigma <- function(chain) {
  cols <- grep("^sigma", colnames(chain), value = TRUE)
  return(chain[, cols])
}

subset_psi <- function(chain) {
  cols <- grep("^psi", colnames(chain), value = TRUE)
  return(chain[, cols])
}

subset_mu <- function(chain) {
  cols <- grep("^mu", colnames(chain), value = TRUE)
  return(chain[, cols])
}

subset_x <- function(chain) {
  cols <- grep("^x", colnames(chain), value = TRUE)
  return(chain[, cols])
}


# CRPS function -----------------------------------------------------------

# CRPS from posterior predictive samples (numeric vector)
# and one observed value.
crps_from_samples <- function(samples,y) {
  samples <- as.numeric(samples)
  M <- length(samples)
  if (M < 2L) stop("Need at least 2 samples")
  
  # Term 1: E|Y - y|  (this is just MAE of the samples vs truth).
  term1 <- mean(abs(samples - y))
  
  # Term 2: E|Y - Y'|, efficient computation using sorted samples.
  y_sorted <- sort(samples)
  coeffs <- (2 * seq_len(M) - M - 1)
  T2 <- (2 / (M^2)) * sum(coeffs * y_sorted)
  
  # CRPS = term1 - 0.5 * term2
  term1 - 0.5 * T2
}
