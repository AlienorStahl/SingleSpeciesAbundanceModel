# Simulations - running the various models ####

# Init ----
  # cleaning the environment
  rm(list=ls())
  # restarting R session
  .rs.restartR()
  
  # Libraries
  library(mvabund)      # for manyglm function
  library(foreach)      # for parallel loop
  library(doSNOW)       # for parallel coding
  library(ecoCopula)    # for gaussian copulas
  
  # starting the parallel coding
  cl <- makeSOCKcluster(4, outfile="test2.txt")
  registerDoSNOW(cl)
  
  # Importing all functions
  source("Scripts/Simulations_functions.R")
  
# end Init
  
# I. Setting the parameters ----
  n.universes = 30 # number of universes/set of parameters
  universes = 1:n.universes
  n.replicates = 10 # number of samples for each universe
  n.species = 20 # number of species in each universe
  n.lakes.landscape.interest <- 1000 # number of sites
  sample.size = 500 # number of sites used to fit the models

  
  
  # Formulas of the models
  # formulas with only the environment
  X_formulas = c(as.formula("Y.fit ~ X1 + X2 + X3")
                 , as.formula("Y.fit ~ X1 + X2")
                 , as.formula("Y.fit ~ X1 + X3")
                 , as.formula("Y.fit ~ X2 + X3")
                 , as.formula("Y.fit ~ X1")
                 , as.formula("Y.fit ~ X2")
                 , as.formula("Y.fit ~ X3")
  )
  # formulas with the latent variables
  L_formulas <- c(as.formula("Y.fit ~ L")
                  , as.formula("Y.fit ~ X1 + X2 + L")
                  , as.formula("Y.fit ~ X1 + X3 + L")
                  , as.formula("Y.fit ~ X2 + X3 + L")
                  , as.formula("Y.fit ~ X1 + L")
                  , as.formula("Y.fit ~ X2 + L")
                  , as.formula("Y.fit ~ X3 + L")
  )
  
  
  # List to record the results
  results.list.original <- rep(list(matrix(NA,n.replicates,n.species)), length(X_formulas)+length(L_formulas))
  list.names <- c("X123", "X12", "X13", "X23", "X1", "X2", "X3"
                  , "L", "X12L", "X13L", "X23L", "X1L", "X2L", "X3L")
  names(results.list.original) <- list.names
  
  
  sensitivity <- specificity <- TSS <- MAPE <- RMSPE <- RMSE <- SMAPE <- RMRPE <- pred.bias <- results.list.original
  oracle.sensitivity <- oracle.specificity <- oracle.TSS <- oracle.MAPE <- oracle.RMSPE <- oracle.RMSE <- oracle.SMAPE <- oracle.RMRPE <- oracle.pred.bias <- matrix(NA,n.replicates,n.species)
  colnames(oracle.sensitivity) <- colnames(oracle.specificity) <- colnames(oracle.TSS) <- colnames(oracle.MAPE) <- colnames(oracle.RMSPE) <- colnames(oracle.RMSE) <- colnames(oracle.SMAPE) <- colnames(oracle.RMRPE) <- colnames(oracle.pred.bias) <- paste0("Spc", 1:n.species)
  
  # lists to record info about abund and occupancy
  occupancy <- abundance_sd <- abundance_tot <- abundance_mean <- matrix(NA, n.replicates, n.species)
  colnames(occupancy) <- colnames(abundance_sd) <- colnames(abundance_tot) <- colnames(abundance_mean) <- paste0("Spc", 1:n.species)
  

  
# end I.
  
# II. General model (sample.size=500)----

  set.seed(5525)
  seeds.universes <- sample(1:10000, n.universes)
  
  mult_uni(universes
           , seeds.universes
           , n.replicates
           , n.species
           , n.lakes.landscape.interest
           , sample.size
           , X_formulas
           , L_formulas)
  
  stopCluster(cl)
  
# end II.
  
# III. All other models----
  # A. 100 sites ----
  set.seed(5525)
  seeds.universes <- sample(1:10000, n.universes)
  
  mult_uni(universes
           , seeds.universes
           , n.replicates
           , n.species
           , n.lakes.landscape.interest
           , 100
           , X_formulas
           , L_formulas)
  
  # B. 200 sites ----
  set.seed(5525)
  seeds.universes <- sample(1:10000, n.universes)
  
  mult_uni(universes
           , seeds.universes
           , n.replicates
           , n.species
           , n.lakes.landscape.interest
           , 200
           , X_formulas
           , L_formulas)
  
  # C. 300 sites ----
  set.seed(5525)
  seeds.universes <- sample(1:10000, n.universes)
  
  mult_uni(universes
           , seeds.universes
           , n.replicates
           , n.species
           , n.lakes.landscape.interest
           , 300
           , X_formulas
           , L_formulas)
  
  stopCluster(cl)
  
# end ----
  
  
  