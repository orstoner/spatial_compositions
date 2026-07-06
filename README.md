# Software for "A Bayesian hierarchical model for spatial compositional data with zeros and missing values", by 
Catherine Holland, Oliver Stoner, and Tereza Neocleous

This is the R code used to produce all the results in the paper,
however we do not have permission to share the original data
and instead provide a simulated compositional dataset
with similar spatial patterns.

The file "simulated_tree_proportions.rds" contains the simulated
tree species data produced by "simulate_data.R". It is not necessary
to run this script but it is included for reproducibility or
modification, as desired.

The script "initial_setup.R" prepares the proportions data for modelling
(e.g. centering and scaling the coordinates). It should be run
first, before running any experiments or "simulate_data.R".

We recommend that the three experiment scrips are run sequentially in
different R sessions, due to the large memory usage. System
memory usage can be reduced by increasing the MCMC thinning rate, though
the highest system memory use will likely occur while running 
the Beta-Binomial GAMs in "experiment_1.R".

Several steps involve parallelisation across multiple CPU cores.
Depending on your system setup, you may need to change the number
of cores used inside "makeCluster" and "registerDoParallel".

The script "experiment_1.R" runs the experiment referred to as
"the main experiment" in the paper.
