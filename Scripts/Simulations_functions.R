# Script containing the functions for the simulations ----

# Generate the communities
generate_community <- function(sample.size # number of sites
                               , bX # environmental coefficients for each species
                               , n.species # number of species
                               ) {
  X <- cbind(rep(1, times=sample.size),matrix(rnorm(sample.size * 3, mean=0, sd=1),sample.size,3))
  colnames(X) <- c("b0", "X1","X2","X3")
  Y <- exp(X %*% bX)
  Y <- rpois(sample.size*n.species,Y)
  Y <- matrix(Y,sample.size,n.species,byrow=FALSE)
  colnames(Y) <- paste0("Sp", 1:n.species) 
  return(list(Y=Y, X=X))
}

# calculate the latent variables from PA
copula_lat <- function(PA.universe # presence-absence matrix of each species in each site
                       , species_X # target species index
                       , n.species # number of species
                       ) {
  eco_PA <- ecoCopula::stackedsdm(PA.universe[,-species_X], ~1, data=PA.universe, family="binomial")
  eco_lvs <- ecoCopula::cord(eco_PA, nlv=3)
  eco_lvs <- eco_lvs$scores
  return(eco_lvs)
}

# fitting the environmental models
modfit_X <- function(Y.fit = Y.fit # abundance matrix used to FIT the model
                     , X.fit # environmental matrix used to FIT the model
                     , X.pred # environmental matrix used for PREDICTIONS
                     , X_formulas = X_formulas # formulas of each model to fit
                     ) {
  data.temp <- data.frame(X.fit)
  data.temp$Y.fit=Y.fit
  predX.list <- rep(list(matrix(0,nrow(X.pred),ncol(Y.fit))), length(X_formulas))
  names(predX.list) <- c("X123", "X12", "X13", "X23", "X1", "X2", "X3")
  for(k in 1:length(X_formulas)) {
    f_temp <- X_formulas[[k]]
    # fit the model on the sample
    mod <- manyglm(formula = as.formula(f_temp), family="poisson", data=data.temp) 
    # predicted values based on the sample model and the new X (rest of universe)
    # this is the conditional expectation
    predX.list[[k]] <- predict(mod, newdata = data.frame(X.pred), type="response")
  }
  return(predX.list)
}

# fitting the latent models
modfit_L <- function(fit # dataframe containing the env, abundance and latent matrix to FIT the model
                     , pred # dataframe containing the env, and latent matrix for PREDICTIONS
                     , L_formulas # formula of each model to fit
                     , predL.list # matrix to record the results
                     , species_X # target species
                     , lvs # matrix of the latent variables
                     ) {
  for (k in 1:length(L_formulas)) {
    f_temp <- gsub(pattern="L", L_formulas[k], replacement=paste(lvs, collapse=" + "))
    mod <- manyglm(formula=as.formula(f_temp),family="poisson", data=fit) 
    predL.list[,k] <- predict(mod, newdata=pred,type="response")
  }
  predL.list[,"Species"] <- species_X
  return(predL.list)
}

# function to calculate various error metrics
abundance.predictive.errors <- function(Yi, # OBSERVED value of abundance
                                        Pi # PREDICTED value of abundance
                                        ){
  # REMINDER: Yi is observed and Pi is predicted
  # checking how many sites contain present species
  n <- sum(!is.na(Yi))
  # error metrics
  MAPE <- (1/n) * sum(abs(Yi - Pi) / Yi, na.rm=TRUE) * 100 # Mean absolute percentage error (MAPE)
  RMSPE <- sqrt((1/n) * sum((Yi - Pi)/Yi, na.rm=TRUE)^2) * 100 # Root mean squared percentage error (RMSPE)
  RMSE <- sqrt(1/n * sum((Yi - Pi)^2 / Yi^2, na.rm=TRUE)) * 100 # Relative mean squared error (RMSE)
  SMAPE <- (1/n) * sum(abs(Yi - Pi) / (abs(Yi) + abs(Pi)), na.rm=TRUE) * 100 # Symmetric mean absolute percentage error (SMAPE)
  RMRPE <- sqrt(1/n * sum(log(Pi/Yi)^2, na.rm=TRUE)) * 100
  pred.bias <- 1/n * sum(log(Pi/Yi), na.rm=TRUE)
  result <- list(MAPE=MAPE,RMSPE=RMSPE,RMSE=RMSE,SMAPE=SMAPE,RMRPE=RMRPE,pred.bias=pred.bias)
  return(result)
}

# Complete simulation function
# generate communities, fit environmental and latent models, calculate error metrics
mult_uni <- function(universes # vector with number of the landscape to run
                     , seeds.universes # vector of seed for each landscape
                     , n.replicates # number of replicates per landscape
                     , n.species # number of species
                     , n.lakes.landscape.interest # total number of sites in the landscape
                     , sample.size # number of sites used to fit the model
                     , X_formulas # formulas for the environmental models
                     , L_formulas # formulas for the latent models
                     ) {
  # universe = landscape
  # loop universe ----
  for (i in universes) { 
    set.seed(seeds.universes[i])
    # setting the parameters of the universe
    bX <- rbind(runif(n.species, -2.4, 1.2), matrix(runif(n.species*3,-0.8,0.8),3,n.species))
    colnames(bX) <- paste0("Spc", 1:n.species)
    rownames(bX) <- c("b0", "b1", "b2", "b3")
    # loop replicates ----
    for (j in 1:n.replicates) { 
      print(paste0(i, ".", j))
      # creating the sites for fitting
      # .fit means it's a dataset used to FIT models
      universe.fit <- generate_community(sample.size, bX, n.species)
      X.fit <- universe.fit$X #environmental variables
      Y.fit <- universe.fit$Y # abundance
      PA.fit <- ifelse(Y.fit >= 1, 1, 0) # transform into presence-absence
      
      ## simulate a new X and generate the remaining samples in the universe,i.e.,
      # n.lakes.landscape.interest - sample.size
      # .pred means it's a dataset used to make PREDICTIONS or compare with PREDICTIONS
      remaining.universe.size <- n.lakes.landscape.interest - sample.size
      universe.pred <- generate_community(remaining.universe.size, bX, n.species)
      X.pred <- universe.pred$X # env
      Y.pred <- universe.pred$Y # abundance
      PA.pred <- ifelse(Y.pred >= 1, 1, 0) # trasnform into presence absence
      # Get PA for whole universe
      PA.universe <- rbind(PA.fit, PA.pred)
      
      # fit the environmental models ----
      # this gives back the conditional expectations
      predX.list <- modfit_X(Y.fit = Y.fit
                             , X.fit = X.fit
                             , X.pred = X.pred
                             , X_formulas = X_formulas)

      # Running the latent models
      # table to save the results
      mat_base <- matrix(NA,nrow(Y.pred),length(L_formulas)+1)
      colnames(mat_base) <- c("L", "X12L", "X13L", "X23L", "X1L", "X2L", "X3L", "Species")
      predL.list <- mat_base
      
      
      # loop species ----
      results <- foreach(species_X = 1:n.species, .packages = c("mvabund", "glmnet", "reshape2", "rlang", "ecoCopula")
                         , .export=c("copula_lat", "modfit_L")) %dopar% {
                           # get the latent variables
                           eco_lvs <- copula_lat(PA.universe, species_X, n.species)
                           # make a table with environment + latent for fitting
                           temp.fit <- data.frame(Y.fit=Y.fit[,species_X]
                                                  , X.fit
                                                  , eco_lvs[1:sample.size,])
                           colnames(temp.fit) <- c("Y.fit", colnames(X.fit), colnames(eco_lvs))
                           # make a table with environment + latent for predictions
                           temp.pred <- data.frame(X.pred
                                                   , eco_lvs[(sample.size+1):n.lakes.landscape.interest,])
                           colnames(temp.pred) <- c(colnames(X.pred), colnames(eco_lvs))
                           # fit latent model on original data + predict
                           # this gives back the conditional expectations
                           list <- modfit_L(temp.fit, temp.pred, L_formulas, predL.list, species_X, colnames(eco_lvs))
                           
                           results <- list
                           return(results)
                         }
      # results formatting ----
      # formatting the results from the latent models
      predL.list <- rep(list(matrix(NA,nrow(Y.pred),n.species)), length(L_formulas))
      names(predL.list) <- list.names[8:14]
      for (k in 1:length(L_formulas)) {
        predL.list[[k]] <- sapply(1:length(results), function(x) results[[x]][,k])
      }
      # putting together the environmental models and the latent models
      predFM.list <- c(predX.list, predL.list)
      
      
      # oracle model (conditional expectations)
      oracle <- matrix((exp(X.pred %*% bX)),remaining.universe.size,n.species,byrow=FALSE)
      colnames(oracle) <- paste0("Sp", 1:n.species)
      # correcting the mean for the positive values
      # this is the equivalent of 1-ppois(0, pred) for the other models
      # aka what's the probability of not having a 0
      oracle.probpresence <- 1-exp(-oracle)
      
      # calculate the threshold per species
      threshold.species <- apply(PA.fit, 2, sum)/sample.size
      
      # empty tables to save the results
      oracle.PA <- matrix(data=NA, nrow=remaining.universe.size, ncol=n.species)
      PA.predicted <- matrix(data=NA, nrow=remaining.universe.size, ncol=n.species)
      
      
      # for lakes that are expected to be present, those are the expected abundance
      # conditional expectations given the lake is not empty
      oracle.pred.abundance.given.present = oracle / oracle.probpresence
      # keep only sites where species are truly present
      # aka remove sites where species are truly absent from the predicted abundance
      oracle.pred.abundance.given.present[PA.pred==0] <- NA
      formated_Ypred <- ifelse(PA.pred==0, NA, Y.pred)
      
      # metrics for the oracle
      for (index in 1:n.species) {
        # get the PA predicted by the oracle
        oracle.PA[,index] <- ifelse(oracle.probpresence[,index] > threshold.species[index], 1, 0)
        
        oracle.results.species <- abundance.predictive.errors(formated_Ypred[,index], oracle.pred.abundance.given.present[,index])
        oracle.MAPE[j,index] <- oracle.results.species$MAPE
        oracle.RMSPE[j,index] <- oracle.results.species$RMSPE
        oracle.RMSE[j,index] <- oracle.results.species$RMSE
        oracle.SMAPE[j,index] <- oracle.results.species$SMAPE
        oracle.RMRPE[j,index] <- oracle.results.species$RMRPE
        oracle.pred.bias[j,index] <- oracle.results.species$pred.bias
      }
      
      # and then calculate PA metrics
      oracle.TP <- apply(oracle.PA==1&PA.pred==1, 2, sum) / apply(PA.pred==1, 2, sum)
      oracle.FP <- apply(oracle.PA==1&PA.pred==0, 2, sum) / apply(PA.pred==0, 2, sum)
      oracle.TN <- apply(oracle.PA==0&PA.pred==0, 2, sum) / apply(PA.pred==0, 2, sum)
      oracle.FN <- apply(oracle.PA==0&PA.pred==1, 2, sum) / apply(PA.pred==1, 2, sum)
      # and now the PA index
      oracle.sensitivity_temp = (oracle.TP / (oracle.TP + oracle.FN))
      oracle.specificity_temp = (oracle.TN / (oracle.FP + oracle.TN)) 
      oracle.TSS_temp = oracle.sensitivity_temp + oracle.specificity_temp -1
      # saving it for the replicate
      oracle.sensitivity[j,] <- oracle.sensitivity_temp
      oracle.specificity[j,] <- oracle.specificity_temp
      oracle.TSS[j, ] <- oracle.TSS_temp
      
      # metrics for the other models
      for (index in list.names) {
        # get the probability that the value predicted is different from 0
        pred.presence = 1-ppois(0, predFM.list[[index]])
        for (species in 1:n.species) {
          # here use the threshold per species to see if each site is present or absent
          # if above the threshold, present, else absent
          PA.predicted[,species] <- ifelse(pred.presence[,species] > threshold.species[species], 1, 0)
        }
        # and then calculate PA metrics
        TP.metric <- apply(PA.predicted==1&PA.pred==1, 2, sum) / apply(PA.pred==1, 2, sum)
        FP.metric <- apply(PA.predicted==1&PA.pred==0, 2, sum) / apply(PA.pred==0, 2, sum)
        TN.metric <- apply(PA.predicted==0&PA.pred==0, 2, sum) / apply(PA.pred==0, 2, sum)
        FN.metric <- apply(PA.predicted==0&PA.pred==1, 2, sum) / apply(PA.pred==1, 2, sum)
        
        # and now the PA index
        sensitivity[[index]][j,] = (TP.metric / (TP.metric + FN.metric))
        specificity[[index]][j,] = (TN.metric / (FP.metric + TN.metric))
        TSS[[index]][j,] = sensitivity[[index]][j,] + specificity[[index]][j,] - 1
        
        
        # for lakes that are expected to be present, those are the expected abundance
        # conditional expectations given the lake is not empty
        pred.abundance.given.present = predFM.list[[index]] / pred.presence
        # keep only sites where species are truly present
        # aka remove sites where species are truly absent from the predicted abundance
        pred.abundance.given.present[PA.pred==0] <- NA
        # removing sites where species is truly absent from our TRUE abundance
        formated_Ypred <- ifelse(PA.pred==0, NA, Y.pred)
        # calculating the various metrics on abundance Pedro mentioned
        for (species in 1:n.species) {
          results.species <- abundance.predictive.errors(formated_Ypred[,species], pred.abundance.given.present[,species])
          MAPE[[index]][j,species] <- results.species$MAPE
          RMSPE[[index]][j,species] <- results.species$RMSPE
          RMSE[[index]][j,species] <- results.species$RMSE
          SMAPE[[index]][j,species] <- results.species$SMAPE
          RMRPE[[index]][j,species] <- results.species$RMRPE
          pred.bias[[index]][j,species] <- results.species$pred.bias
        }
      }
      
      
      
      # calculating the occurance and abundance of each species
      occupancy[j,] <- apply(PA.universe, 2, sum)
      abundance_tot[j,] <- apply(Y.pred, 2, sum)+ apply(Y.fit, 2, sum)
      abundance_sd[j,] <- apply(rbind(Y.pred, Y.fit), 2, sd)
    }
    # getting a file name where we'll save the results
    vector_name <- paste0("universe_", i)
    file_name <- paste0("Results/Sample.size_", sample.size, "/", vector_name, ".RData")
    # saving the important results
    save(bX, occupancy, abundance_sd, abundance_tot
         , sensitivity, specificity, oracle.sensitivity, oracle.specificity
         , TSS
         , MAPE, RMSPE, RMSE, SMAPE, RMRPE, pred.bias, oracle.RMSE
         , oracle.TSS, oracle.MAPE, oracle.RMSPE, oracle.RMSE, oracle.SMAPE
         , oracle.RMRPE, oracle.pred.bias
         , file = file_name)
    
  }
  # closing the cluster in case it's left to run on its own to avoid rogue sessions
  stopCluster(cl)
  
}
