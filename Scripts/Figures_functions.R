# Script containing the functions for the plots ----

# General functions ----

# function to generate the communities
generate_community <- function(sample.size # number of sites in the landscape
                               , X # matrix of environmental values
                               , bX # coefficients linking species and env variables
                               , n.species # number of species in the landscape
) {
  Y <- exp(X %*% bX)
  Y <- rpois(sample.size*n.species,Y)
  Y <- matrix(Y,sample.size,n.species,byrow=FALSE)
  colnames(Y) <- paste0("Sp", 1:n.species) 
  return(list(Y=Y, X=X))
}

# function to calculate the gaussian copulas
copula_lat <- function(PA.universe # matrix containing the presence absence for each species in each sites
) {
  eco_PA <- ecoCopula::stackedsdm(PA.universe, ~1, data=PA.universe, family="binomial")
  eco_lvs <- ecoCopula::cord(eco_PA, nlv=1)
  eco_lvs <- eco_lvs$scores
  return(eco_lvs)
}

# function to calculate the various error metrics picked
abundance.predictive.errors <- function(Yi # OBSERVED value of abundance
                                        ,Pi # PREDICTED value of abundance
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
  # formatting results into a list
  result <- list(MAPE=MAPE,RMSPE=RMSPE,RMSE=RMSE,SMAPE=SMAPE,RMRPE=RMRPE,pred.bias=pred.bias)
  return(result)
}

# Formatting functions for plots ----

# function to retrieve the coefficients for each species and assign which model is considered as which
# to use in a loop (for each universe)
bX_transformation <- function(universe.number # universe number considered
                              , n.species # number of species
                              , n.sites # number of sites used to fit the data
                              ) {
  file_name <- paste0("Results/Sample.size_", n.sites, "/Simul_universe_", universe.number, ".RData")
  load(file_name)
  # seeing which parameter is the most important (higher value)
  bX_ordered <- data.frame(t(apply(-abs(bX[c("b1", "b2", "b3"),]), 2, rank)))
  colnames(bX_ordered) <- c("b1", "b2", "b3")
  # Setting most important parameter as "high" performing model with a single variable
  bX_summary <- data.frame(ifelse(bX_ordered == 1, "Mod.1V.high", ifelse(bX_ordered==2, "Mod.1V.mid", "Mod.1V.low")))
  bX_summary$Species <- paste0("Spc", 1:n.species)
  colnames(bX_summary) <- c("X1", "X2", "X3", "Species")
  # combining high and mid performing models with 1 var into a high performing model with 2 variables
  # high and low performing models with 1 var into a mid performing model with 2 variables
  # mid and low performing models with 1 var into a low performing model with 2 variables
  bX_summary$X12 <- ifelse((grepl("high", paste0(bX_summary$X1, bX_summary$X2))+grepl("mid", paste0(bX_summary$X1, bX_summary$X2))==2)==TRUE
                           , "Mod.2V.high"
                           , ifelse((grepl("high", paste0(bX_summary$X1, bX_summary$X2))+grepl("low", paste0(bX_summary$X1, bX_summary$X2))==2)==TRUE
                                    , "Mod.2V.mid"
                                    , "Mod.2V.low"))
  bX_summary$X13 <- ifelse((grepl("high", paste0(bX_summary$X1, bX_summary$X3))+grepl("mid", paste0(bX_summary$X1, bX_summary$X3))==2)==TRUE
                           , "Mod.2V.high"
                           , ifelse((grepl("high", paste0(bX_summary$X1, bX_summary$X3))+grepl("low", paste0(bX_summary$X1, bX_summary$X3))==2)==TRUE
                                    , "Mod.2V.mid"
                                    , "Mod.2V.low"))
  bX_summary$X23 <- ifelse((grepl("high", paste0(bX_summary$X3, bX_summary$X2))+grepl("mid", paste0(bX_summary$X3, bX_summary$X2))==2)==TRUE
                           , "Mod.2V.high"
                           , ifelse((grepl("high", paste0(bX_summary$X3, bX_summary$X2))+grepl("low", paste0(bX_summary$X3, bX_summary$X2))==2)==TRUE
                                    , "Mod.2V.mid"
                                    , "Mod.2V.low"))
  bX_summary$X123 <- "Benchmark"
  # Adding the info for models with the latence
  bX_latent <- data.frame(sapply(colnames(bX_summary)[colnames(bX_summary)!="Species"], function(x) paste0(bX_summary[,x], "L")))
  colnames(bX_latent) <- paste0(colnames(bX_summary)[colnames(bX_summary)!="Species"], "L")
  bX_latent$Species <- paste0("Spc", 1:n.species)
  bX_sumall <- dplyr::left_join(bX_summary, bX_latent, by="Species")
  bX_sumall$L <- "Mod.L"
  # melting it into a simple frame we can then join to the rest of the data
  bX_sumall <- melt(bX_sumall, id.vars="Species")
  colnames(bX_sumall) <- c("Species", "Model", "GenModel")
  return(bX_sumall)
}

# function to retrieve and format the metrics for a certain sample size
format_metrics <- function(n.species # number of species
                         , universes # vector of universes to consider
                         , n.sites # number of sites used to fit the data
                         ) {
  data_all <- data.frame(matrix(data=NA, nrow = 0, ncol = 16))
  colnames(data_all) <- c("Universe", "Species", "Model", "GenModel"
                          , "TSS", "SEDI", "MAPE", "RMSPE"
                          , "RMSE", "SMAPE", "RMRPE", "pred.bias", "RMSE_vOracle"
                          , "Avg_Abund.tot", "Avg_Abund.mean", "Avg_occupancy")
  for (i in universes) {
    # loading all data
    file_name <- paste0("Results/Sample.size_", n.sites, "/Simul_universe_", i, ".RData")
    load(file_name)
    
    # formating the data on abundance and PA
    all_data <- data.frame(Species = paste0("Spc", 1:n.species)
                           , Avg_Abund.tot = NA
                           , Avg_Abund.mean = NA
                           , Avg_occupancy = NA)
    all_data$Avg_Abund.tot <- apply(abundance_tot, 2, mean)
    all_data$Avg_occupancy <- apply(occupancy, 2, mean)
    all_data$Avg_Abund.mean <- apply(abundance_tot, 2, mean)/1000
    
    # formatting the metrics
    temp_TSS <- sapply(1:14, function(x) apply(TSS[[x]], 2, mean))
    temp_sensitivity <- sapply(1:14, function(x) apply(sensitivity[[x]], 2, mean))
    temp_specificity <- sapply(1:14, function(x) apply(specificity[[x]], 2, mean))
    # temp_SEDI <- sapply(1:14, function(x) apply(SEDI[[x]], 2, mean))
    temp_MAPE <- sapply(1:14, function(x) apply(MAPE[[x]], 2, mean))
    temp_RMSPE <- sapply(1:14, function(x) apply(RMSPE[[x]], 2, mean))
    temp_RMSE <- sapply(1:14, function(x) apply(RMSE[[x]], 2, mean))
    temp_SMAPE <- sapply(1:14, function(x) apply(SMAPE[[x]], 2, mean))
    temp_RMRPE <- sapply(1:14, function(x) apply(RMRPE[[x]], 2, mean))
    temp_pred.bias <- sapply(1:14, function(x) apply(pred.bias[[x]], 2, mean))
    
    temp_TSS_vOracle <- sapply(1:14, function(x) apply(TSS[[x]]/oracle.TSS, 2, mean))
    temp_sensitivity_vOracle <- sapply(1:14, function(x) apply(sensitivity[[x]]/oracle.sensitivity, 2, mean))
    temp_specificity_vOracle <- sapply(1:14, function(x) apply(specificity[[x]]/oracle.specificity, 2, mean))
    temp_MAPE_vOracle <- sapply(1:14, function(x) apply(MAPE[[x]]/oracle.MAPE, 2, mean))
    temp_RMSPE_vOracle <- sapply(1:14, function(x) apply(RMSPE[[x]]/oracle.RMSPE, 2, mean))
    temp_RMSE_vOracle <- sapply(1:14, function(x) apply(RMSE[[x]]/oracle.RMSE, 2, mean))
    temp_SMAPE_vOracle <- sapply(1:14, function(x) apply(SMAPE[[x]]/oracle.SMAPE, 2, mean))
    temp_RMRPE_vOracle <- sapply(1:14, function(x) apply(RMRPE[[x]]/oracle.RMRPE, 2, mean))
    temp_pred.bias_vOracle <- sapply(1:14, function(x) apply(pred.bias[[x]] - oracle.pred.bias, 2, mean))
    
    
    colnames(temp_TSS) <- colnames(temp_sensitivity) <- colnames(temp_specificity) <- colnames(temp_MAPE) <- colnames(temp_RMSPE) <- colnames(temp_RMSE) <- colnames(temp_SMAPE) <- colnames(temp_RMRPE) <- colnames(temp_pred.bias) <- names(TSS)[1:14]
    colnames(temp_TSS_vOracle) <- colnames(temp_sensitivity_vOracle) <- colnames(temp_specificity_vOracle) <- colnames(temp_MAPE_vOracle) <- colnames(temp_RMSPE_vOracle) <- colnames(temp_RMSE_vOracle) <- colnames(temp_SMAPE_vOracle) <- colnames(temp_RMRPE_vOracle) <- colnames(temp_pred.bias_vOracle) <- names(TSS)[1:14]
    rownames(temp_TSS) <- rownames(temp_sensitivity) <- rownames(temp_specificity) <- rownames(temp_MAPE) <- rownames(temp_RMSPE) <- rownames(temp_RMSE) <- rownames(temp_SMAPE) <- rownames(temp_RMRPE) <- rownames(temp_pred.bias) <- paste0("Spc", 1:n.species)
    rownames(temp_TSS_vOracle) <- rownames(temp_sensitivity_vOracle) <- rownames(temp_specificity_vOracle) <- rownames(temp_MAPE_vOracle) <- rownames(temp_RMSPE_vOracle) <- rownames(temp_RMSE_vOracle) <- rownames(temp_SMAPE_vOracle) <- rownames(temp_RMRPE_vOracle) <- rownames(temp_pred.bias_vOracle) <- paste0("Spc", 1:n.species)
    
    list_temp <-list(TSS = temp_TSS
                     , sensitivity = temp_sensitivity
                     , specificity = temp_sensitivity
                     , MAPE = temp_MAPE
                     , RMSPE = temp_RMSPE
                     , RMSE = temp_RMSE
                     , SMAPE = temp_SMAPE
                     , RMRPE = temp_RMRPE
                     , pred.bias = temp_pred.bias
                     , TSS_vOracle = temp_TSS_vOracle
                     , sensitivity_vOracle = temp_sensitivity_vOracle
                     , specificity_vOracle = temp_specificity_vOracle
                     , MAPE_vOracle = temp_MAPE_vOracle
                     , RMSPE_vOracle = temp_RMSPE_vOracle
                     , RMSE_vOracle = temp_RMSE_vOracle
                     , SMAPE_vOracle = temp_SMAPE_vOracle
                     , RMRPE_vOracle = temp_RMRPE_vOracle
                     , pred.bias_vOracle = temp_pred.bias_vOracle
                     )
    names_temp <- c("TSS", "sensitivity", "specificity", "MAPE", "RMSPE", "RMSE", "SMAPE", "RMRPE", "pred.bias"
                    , "TSS_vOracle", "sensitivity_vOracle", "specificity_vOracle", "MAPE_vOracle", "RMSPE_vOracle", "RMSE_vOracle", "SMAPE_vOracle", "RMRPE_vOracle", "pred.bias_vOracle")
    test <- lapply(names_temp, function(x) melt(list_temp[[x]], value.name=x, varnames=c("Species", "Model")))
    names(test) <- names_temp
    test <- Reduce(merge, test)
    
    # formating the new names of the models
    bX_all <- bX_transformation(i, n.species, n.sites)
    
    # joining all the information
    test <- dplyr::left_join(test, bX_all, by=c("Species", "Model"))
    data_uni_temp <- test
    data_uni_temp <- dplyr::left_join(test, all_data, by="Species")
    data_uni_temp$Universe <- i
    data_all <- rbind(data_all, data_uni_temp)
  }
  return(data_all)
}

# function to format the predictions of the oracle model
format_oracle <- function(n.species # number of species
                          , universes # vector with which universe to consider
                          , n.sites # number of sites used to fit the data
                          ) {
  data_all <- data.frame(matrix(data=NA, nrow = 0, ncol = 16))
  colnames(data_all) <- c("Universe", "Species"
                          , "TSS", "sensitivity", "specificity", "MAPE", "RMSPE"
                          , "RMSE", "SMAPE", "RMRPE", "pred.bias"
                          , "Avg_Abund.tot", "Avg_Abund.mean", "Avg_occupancy")
  for (i in universes) {
    # loading all data
    file_name <- paste0("Results/Sample.size_", n.sites, "/Simul_universe_", i, ".RData")
    load(file_name)
    
    # formating the data on abundance and PA
    all_data <- data.frame(Species = paste0("Spc", 1:n.species)
                           , Avg_Abund.tot = NA
                           , Avg_Abund.mean = NA
                           , Avg_occupancy = NA)
    all_data$Avg_Abund.tot <- apply(abundance_tot, 2, mean)
    all_data$Avg_occupancy <- apply(occupancy, 2, mean)
    all_data$Avg_Abund.mean <- apply(abundance_sd, 2, mean)
    
    # formatting the metrics
    temp_oracle <- data.frame(TSS = apply(oracle.TSS, 2, mean)
                              , sensitivity = apply(oracle.sensitivity, 2, mean)
                              , specificity = apply(oracle.specificity, 2, mean)
                              , MAPE = apply(oracle.MAPE, 2, mean)
                              , RMSPE = apply(oracle.RMSPE, 2, mean)
                              , RMSE = apply(oracle.RMSE, 2, mean)
                              , SMAPE = apply(oracle.SMAPE, 2, mean)
                              , RMRPE = apply(oracle.RMRPE, 2, mean)
                              , pred.bias = apply(oracle.pred.bias, 2, mean))
    temp_oracle$Species <- paste0("Spc", 1:n.species)
    
    # joining all the information
    data_uni_temp <- dplyr::left_join(temp_oracle, all_data, by="Species")
    data_uni_temp$Universe <- i
    data_all <- rbind(data_all, data_uni_temp)
  }
  data_all$GenModel = "Oracle"
  return(data_all)
}

# function to retrieve and format the metrics for ALL sample size
format_evol <- function(n.species # number of species
                        , universes # vector with which universe to consider
                        ) {
  data_all <- data.frame(matrix(data=NA, nrow = 0, ncol = 17))
  colnames(data_all) <- c("Universe", "Species", "n_sites", "Model", "GenModel"
                          , "TSS", "SEDI", "MAPE", "RMSPE"
                          , "RMSE", "SMAPE", "RMRPE", "pred.bias", "RMSE_vOracle"
                          , "Avg_Abund.tot", "Avg_Abund.mean", "Avg_occupancy")
  files_oi <- dir("Results/")[grep("Sample.size_", dir("Results/"))]
  for (j in files_oi) {
    for (i in universes) {
      # loading all data
      file_name <- paste0("Results/", j, "/Simul_universe_", i, ".RData")
      load(file_name)
      
      # formating the data on abundance and PA
      all_data <- data.frame(Species = paste0("Spc", 1:n.species)
                             , Avg_Abund.tot = NA
                             , Avg_Abund.mean = NA
                             , Avg_occupancy = NA)
      all_data$Avg_Abund.tot <- apply(abundance_tot, 2, mean)
      all_data$Avg_occupancy <- apply(occupancy, 2, mean)
      all_data$Avg_Abund.mean <- apply(abundance_sd, 2, mean)
      
      # formatting the metrics
      temp_TSS <- sapply(1:14, function(x) apply(TSS[[x]], 2, mean))
      temp_sensitivity <- sapply(1:14, function(x) apply(sensitivity[[x]], 2, mean))
      temp_specificity <- sapply(1:14, function(x) apply(specificity[[x]], 2, mean))
      # temp_SEDI <- sapply(1:14, function(x) apply(SEDI[[x]], 2, mean))
      temp_MAPE <- sapply(1:14, function(x) apply(MAPE[[x]], 2, mean))
      temp_RMSPE <- sapply(1:14, function(x) apply(RMSPE[[x]], 2, mean))
      temp_RMSE <- sapply(1:14, function(x) apply(RMSE[[x]], 2, mean))
      temp_SMAPE <- sapply(1:14, function(x) apply(SMAPE[[x]], 2, mean))
      temp_RMRPE <- sapply(1:14, function(x) apply(RMRPE[[x]], 2, mean))
      temp_pred.bias <- sapply(1:14, function(x) apply(pred.bias[[x]], 2, mean))
      
      temp_TSS_vOracle <- sapply(1:14, function(x) apply(TSS[[x]]/oracle.TSS, 2, mean))
      temp_sensitivity_vOracle <- sapply(1:14, function(x) apply(sensitivity[[x]]/oracle.sensitivity, 2, mean))
      temp_specificity_vOracle <- sapply(1:14, function(x) apply(specificity[[x]]/oracle.specificity, 2, mean))
      temp_MAPE_vOracle <- sapply(1:14, function(x) apply(MAPE[[x]]/oracle.MAPE, 2, mean))
      temp_RMSPE_vOracle <- sapply(1:14, function(x) apply(RMSPE[[x]]/oracle.RMSPE, 2, mean))
      temp_RMSE_vOracle <- sapply(1:14, function(x) apply(RMSE[[x]]/oracle.RMSE, 2, mean))
      temp_SMAPE_vOracle <- sapply(1:14, function(x) apply(SMAPE[[x]]/oracle.SMAPE, 2, mean))
      temp_RMRPE_vOracle <- sapply(1:14, function(x) apply(RMRPE[[x]]/oracle.RMRPE, 2, mean))
      temp_pred.bias_vOracle <- sapply(1:14, function(x) apply(pred.bias[[x]] - oracle.pred.bias, 2, mean))
      
      
      colnames(temp_TSS) <- colnames(temp_sensitivity) <- colnames(temp_specificity) <- colnames(temp_MAPE) <- colnames(temp_RMSPE) <- colnames(temp_RMSE) <- colnames(temp_SMAPE) <- colnames(temp_RMRPE) <- colnames(temp_pred.bias) <- names(TSS)[1:14]
      colnames(temp_TSS_vOracle) <- colnames(temp_sensitivity_vOracle) <- colnames(temp_specificity_vOracle) <- colnames(temp_MAPE_vOracle) <- colnames(temp_RMSPE_vOracle) <- colnames(temp_RMSE_vOracle) <- colnames(temp_SMAPE_vOracle) <- colnames(temp_RMRPE_vOracle) <- colnames(temp_pred.bias_vOracle) <- names(TSS)[1:14]
      rownames(temp_TSS) <- rownames(temp_sensitivity) <- rownames(temp_specificity) <- rownames(temp_MAPE) <- rownames(temp_RMSPE) <- rownames(temp_RMSE) <- rownames(temp_SMAPE) <- rownames(temp_RMRPE) <- rownames(temp_pred.bias) <- paste0("Spc", 1:n.species)
      rownames(temp_TSS_vOracle) <- rownames(temp_sensitivity_vOracle) <- rownames(temp_specificity_vOracle) <- rownames(temp_MAPE_vOracle) <- rownames(temp_RMSPE_vOracle) <- rownames(temp_RMSE_vOracle) <- rownames(temp_SMAPE_vOracle) <- rownames(temp_RMRPE_vOracle) <- rownames(temp_pred.bias_vOracle) <- paste0("Spc", 1:n.species)
      
      list_temp <-list(TSS = temp_TSS
                       , sensitivity = temp_sensitivity
                       , specificity = temp_sensitivity
                       , MAPE = temp_MAPE
                       , RMSPE = temp_RMSPE
                       , RMSE = temp_RMSE
                       , SMAPE = temp_SMAPE
                       , RMRPE = temp_RMRPE
                       , pred.bias = temp_pred.bias
                       , TSS_vOracle = temp_TSS_vOracle
                       , sensitivity_vOracle = temp_sensitivity_vOracle
                       , specificity_vOracle = temp_specificity_vOracle
                       , MAPE_vOracle = temp_MAPE_vOracle
                       , RMSPE_vOracle = temp_RMSPE_vOracle
                       , RMSE_vOracle = temp_RMSE_vOracle
                       , SMAPE_vOracle = temp_SMAPE_vOracle
                       , RMRPE_vOracle = temp_RMRPE_vOracle
                       , pred.bias_vOracle = temp_pred.bias_vOracle
      )
      names_temp <- c("TSS", "sensitivity", "specificity", "MAPE", "RMSPE", "RMSE", "SMAPE", "RMRPE", "pred.bias"
                      , "TSS_vOracle", "sensitivity_vOracle", "specificity_vOracle", "MAPE_vOracle", "RMSPE_vOracle", "RMSE_vOracle", "SMAPE_vOracle", "RMRPE_vOracle", "pred.bias_vOracle")
      test <- lapply(names_temp, function(x) melt(list_temp[[x]], value.name=x, varnames=c("Species", "Model")))
      names(test) <- names_temp
      test <- Reduce(merge, test)
      
      # formating the new names of the models
      bX_all <- bX_transformation(i, n.species, 200)
      
      # joining all the information
      test <- dplyr::left_join(test, bX_all, by=c("Species", "Model"))
      data_uni_temp <- test
      data_uni_temp <- dplyr::left_join(test, all_data, by="Species")
      data_uni_temp$Universe <- i
      data_uni_temp$n_sites <- gsub("Sample.size_", "",j)
      data_all <- rbind(data_all, data_uni_temp)
    }
  }
  return(data_all)
}

# Plots for PA ----

# function to plot the correlation of TSS between different models 
plot_PA_correlation <- function(metrics_data # formatted data from format_metrics
                                ) {
  bins_temp <- quantile(metrics_data$Avg_Abund.mean, probs = seq(0, 1, 0.05), digits=7)
  metrics_data$bins_abund <- cut(metrics_data$Avg_Abund.mean, breaks=bins_temp, include.lowest = TRUE)
  bins_temp <- quantile(metrics_data$Avg_occupancy, probs = seq(0, 1, 0.05), digits=7)
  metrics_data$bins_occur <- cut(metrics_data$Avg_occupancy, breaks=bins_temp, include.lowest = TRUE)
  
  data_oracle <- format_oracle(n.species, universes, n.sites)
  bins_temp <- quantile(data_oracle$Avg_Abund.mean, probs = seq(0, 1, 0.05), digits=7)
  data_oracle$bins_abund<- cut(data_oracle$Avg_Abund.mean, breaks=bins_temp, include.lowest = TRUE)
  bins_temp <- quantile(data_oracle$Avg_occupancy, probs = seq(0, 1, 0.05), digits=7)
  data_oracle$bins_occur <- cut(data_oracle$Avg_occupancy, breaks=bins_temp, include.lowest = TRUE)
  
  metrics_data <- metrics_data[metrics_data$GenModel %in% c("Benchmark", "Mod.L"), ]
  metrics_data <- metrics_data[,colnames(metrics_data) %in% colnames(data_oracle)]
  all_data <- rbind(metrics_data, data_oracle)
  bins_df_abund <- data.frame(percentage_abund=seq(5, 100, by=5), bins_abund=levels(all_data$bins_abund))
  bins_df_occur <- data.frame(percentage_occur=seq(5, 100, by=5), bins_occur=levels(all_data$bins_occur))
  all_data <- dplyr::left_join(all_data, bins_df_abund, by="bins_abund")
  all_data <- dplyr::left_join(all_data, bins_df_occur, by="bins_occur")
  all_data$Model <- factor(ifelse(all_data$GenModel=="Mod.L", "Latent", all_data$GenModel)
                           , levels = c("Oracle", "Benchmark", "Latent"))
  
  data_PA <- data.frame(Variable = rep(c("TSS", "Sensitivity", "Specificity"), each = length(all_data$TSS))
                        , Value = c(all_data$TSS, all_data$sensitivity, all_data$specificity)
                        , percentage_occur = c(all_data$percentage_occur, all_data$percentage_occur, all_data$percentage_occur)
                        , Model = c(all_data$Model, all_data$Model, all_data$Model))
  data_PA <- data_PA %>%
    group_by(Variable, percentage_occur, Model) %>%
    summarise(mean_value = mean(Value),
              sd_value = sd(Value,na.rm= TRUE))%>%
    mutate(lower = mean_value - 2*sd_value,
           upper = mean_value + 2*sd_value)
  data_PA$Variable <- factor(data_PA$Variable, levels=c("TSS", "Sensitivity", "Specificity"))
  TSS <- ggplot(data_PA, aes(x=percentage_occur, y=mean_value, group=Model, color=Model))+
    geom_line() +
    geom_pointrange(data = data_PA, aes(ymin=lower, ymax=upper, x=percentage_occur, y=mean_value, color=Model))+
    facet_grid(Variable~., scales = "free_y")+
    theme_bw() +
    theme(panel.background = element_blank(),legend.position = "bottom") +
    scale_color_brewer(palette = "Set1")+
    xlab("Species occurrence percentiles") + ylab("Average value") 
  gridExtra::grid.arrange(TSS)
  ggsave(filename = "Figures/Fig5.jpeg", width = 7.5, height=6, dpi=300, units="in")
  
  
  
}

# plotting the ratio TSS and delta TSS for a fixed number of sites
plot_PA_metrics <- function(metrics_data # formatted data from format_metrics
                            ) {
  # formatting for the plots
  bins_temp <- quantile(metrics_data$Avg_occupancy, probs = seq(0, 1, 0.05), digits=7)
  metrics_data$bins_occur <- cut(metrics_data$Avg_occupancy, breaks=bins_temp, include.lowest = TRUE)
  bins_temp <- quantile(metrics_data$Avg_Abund.mean, probs = seq(0, 1, 0.05), digits=7)
  metrics_data$bins_abund <- cut(metrics_data$Avg_Abund.mean, breaks=bins_temp, include.lowest = TRUE)
  metricsP_data <- metrics_data %>%
    dplyr::group_by(GenModel, bins_occur, bins_abund) %>%
    dplyr::summarise(TSS_vOracle = mean(TSS_vOracle)
                     , sensitivity_vOracle = mean(sensitivity_vOracle)
                     , specificity_vOracle = mean(specificity_vOracle)
                     , MAPE_vOracle = mean(MAPE_vOracle)
                     , RMSPE_vOracle = mean(RMSPE_vOracle)
                     , RMSE_vOracle = mean(RMSE_vOracle)
                     , SMAPE_vOracle = mean(SMAPE_vOracle)
                     , RMRPE_vOracle = mean(RMRPE_vOracle)
                     , pred.bias_vOracle = mean(pred.bias_vOracle))
  bins_df_occur <- data.frame(percentage_occur=seq(5, 100, by=5), bins_occur=levels(metricsP_data$bins_occur))
  bins_df_occur$percentage_occur <- units::set_units(bins_df_occur$percentage_occur, "%")
  metricsP_data <- dplyr::left_join(metricsP_data, bins_df_occur, by="bins_occur")
  bins_df_abund <- data.frame(percentage_abund=seq(5, 100, by=5), bins_abund=levels(metricsP_data$bins_abund))
  bins_df_abund$percentage_abund <- units::set_units(bins_df_abund$percentage_abund, "%")
  metricsP_data <- dplyr::left_join(metricsP_data, bins_df_abund, by="bins_abund")
  library(units)
  metricsP_data$Latent <- "Environment"
  metricsP_data[grep("L", metricsP_data$GenModel),"Latent"] <- "Latent"
  metricsP_data$Latent <- factor(metricsP_data$Latent, levels = c("Environment", "Latent"))
  metricsP_data$DumbModel <- gsub("L", "", metricsP_data$GenModel)
  metricsP_data$DumbModel[metricsP_data$DumbModel=="Mod."] <- "Latent"
  metricsP_data$DumbModel <- gsub("Mod.", "", metricsP_data$DumbModel)
  metricsP_data$n.var <- factor(ifelse(grepl(2, metricsP_data$DumbModel)==TRUE, 2,
                                       ifelse(grepl(1, metricsP_data$DumbModel)==TRUE,1,
                                              ifelse(grepl("Benchmark", metricsP_data$DumbModel)==TRUE, 3, 0)))
                                , levels = c(3,2,1,0))
  metricsP_data$DumbModel <- factor(metricsP_data$DumbModel, levels=c("Latent"
                                                                      , "1V.low", "1V.mid", "1V.high"
                                                                      , "2V.low", "2V.mid", "2V.high", "Benchmark"))
  
  data <- data.frame(Latent = c("Environment", "Latent")
                     , Intercept = c(3.5, 3.5))
  
  metricsP_data <- metricsP_data[!metricsP_data$DumbModel %in% c("Latent", "Benchmark"),]
  
  
  
  # Plotting TSS
  TSS <- ggplot(metricsP_data, aes(percentage_occur, DumbModel)) +
    geom_tile(aes(fill=TSS_vOracle)) +
    theme_bw() +
    theme(panel.background = element_blank()) +
    xlab("Species occurrence percentiles") + ylab("Model") +
    labs(fill='Ratio TSS') +
    scale_fill_viridis_c(limits = c(0, 1), oob=scales::squish, 
                         breaks = c(0,.25,.5,.75,1), 
                         labels = c("<=0","0.25", "0.5", "0.75", ">=1")) +
    facet_grid(Latent~., scales = "free") +
    geom_hline(data=data,aes(yintercept = Intercept))
  
  
  bins_temp <- quantile(metrics_data$Avg_occupancy, probs = seq(0, 1, 0.05), digits=7)
  metrics_data$bins_occur <- cut(metrics_data$Avg_occupancy, breaks=bins_temp, include.lowest = TRUE)
  bins_temp <- quantile(metrics_data$Avg_Abund.mean, probs = seq(0, 1, 0.05), digits=7)
  metrics_data$bins_abund <- cut(metrics_data$Avg_Abund.mean, breaks=bins_temp, include.lowest = TRUE)
  # new column on whether the model has latence or not
  metrics_data$type <- NA
  metrics_data$type[-grep("L", metrics_data$GenModel)] <- "env"
  metrics_data$type[is.na(metrics_data$type)] <- "latent"
  # new column on a general model so we can contrast based on the previous column
  metrics_data$DumbModel <- gsub("L", "", metrics_data$GenModel)
  metrics_data$DumbModel[metrics_data$DumbModel=="Mod."] <- "Benchmark"
  metrics_data <- metrics_data[metrics_data$DumbModel!="Benchmark",]
  compar_data <- metrics_data %>% 
    dplyr::group_by(Universe, Species, DumbModel) %>%
    summarise(delta_TSS = mean(TSS[type=="env"] - TSS[type=="latent"])
              , delta_sensitivity = mean(sensitivity[type=="env"] - sensitivity[type=="latent"])
              , delta_specificity = mean(specificity[type=="env"] - specificity[type=="latent"])
              , delta_MAPE = mean(MAPE[type=="env"] - MAPE[type=="latent"])
              , delta_RMSPE = mean(RMSPE[type=="env"] - RMSPE[type=="latent"])
              , delta_RMSE = mean(RMSE[type=="env"] - RMSE[type=="latent"])
              , delta_SMAPE = mean(SMAPE[type=="env"] - SMAPE[type=="latent"])
              , delta_RMRPE = mean(RMRPE[type=="env"] - RMRPE[type=="latent"])
              , delta_pred.bias = mean(abs(pred.bias[type=="env"]) - abs(pred.bias[type=="latent"]))
              , Avg_occupancy = unique(Avg_occupancy)
              , bins_occur = unique(bins_occur)
              , bins_abund = unique(bins_abund))
  comparP_data <- compar_data %>%
    dplyr::group_by(DumbModel, bins_abund, bins_occur) %>%
    dplyr::summarise(delta_TSS = mean(delta_TSS)
                     , delta_sensitivity = mean(delta_sensitivity)
                     , delta_specificity = mean(delta_specificity)
                     , delta_MAPE = mean(delta_MAPE)
                     , delta_RMSPE = mean(delta_RMSPE)
                     , delta_RMSE = mean(delta_RMSE)
                     , delta_SMAPE = mean(delta_SMAPE)
                     , delta_RMRPE = mean(delta_RMRPE)
                     , delta_pred.bias = mean(delta_pred.bias))
  bins_df_occur <- data.frame(percentage_occur=seq(5, 100, by=5), bins_occur=levels(comparP_data$bins_occur))
  bins_df_occur$percentage_occur <- units::set_units(bins_df_occur$percentage_occur, "%")
  comparP_data <- dplyr::left_join(comparP_data, bins_df_occur, by="bins_occur")
  bins_df_abund <- data.frame(percentage_abund=seq(5, 100, by=5), bins_abund=levels(comparP_data$bins_abund))
  bins_df_abund$percentage <- units::set_units(bins_df_abund$percentage_abund, "%")
  comparP_data <- dplyr::left_join(comparP_data, bins_df_abund, by="bins_abund")
  
  comparP_data$DumbModel <- gsub("Mod.", "", comparP_data$DumbModel)
  comparP_data$DumbModel <- factor(comparP_data$DumbModel, levels=c("1V.low", "1V.mid", "1V.high"
                                                                    , "2V.low", "2V.mid", "2V.high"))
  
  
  TSS_delta <- ggplot(comparP_data, aes(percentage_occur, DumbModel)) +
    geom_tile(aes(fill=delta_TSS)) +
    theme_bw()+
    theme(panel.background = element_blank()) +
    xlab("Species occurrence percentiles") + ylab("Model") +
    labs(fill='\u0394 TSS') +    
    scale_fill_gradient2(high = "blue", mid="white", low = "red", oob=scales::squish) +
    geom_hline(data=data,aes(yintercept = 3.5))
  
  cowplot::plot_grid(TSS, TSS_delta, nrow=2, align = "v", axis = "Model", rel_heights = c(1.5,1))
  ggsave(filename = "Figures/Fig4.jpeg", width = 8, height=9, dpi=600, units="in")
}

# plotting the ratio TSS, sensitivity and specificity depending on number of species, sites and average abundance
plot_PA_sites <- function(metrics_evol # formatted data from format_evol
                         ) {
  metrics_evol <- metrics_evol[metrics_evol$GenModel %in% c("Benchmark", "Mod.L", "Mod.1V.high", "Mod.2V.high"),]
  metrics_evol <- metrics_evol[,colnames(metrics_evol) %in% c("Species", "TSS_vOracle", "sensitivity_vOracle", 
                                                              "specificity_vOracle", "GenModel"
                                                              , "Avg_occupancy", "Universe", "n_sites")]
  bins_temp <- quantile(metrics_evol$Avg_occupancy, probs = seq(0, 1, 0.05), digits=7)
  metrics_evol$bins_occur <- cut(metrics_evol$Avg_occupancy, breaks=bins_temp, include.lowest = TRUE)
  bins_df_occur <- data.frame(percentage_occur=seq(5, 100, by=5), bins_occur=levels(metrics_evol$bins_occur))
  metrics_evol <- dplyr::left_join(metrics_evol, bins_df_occur, by="bins_occur")
  
  data_all <- data.frame(Universe = c(metrics_evol$Universe, metrics_evol$Universe, metrics_evol$Universe)
                         , Avg_occupancy = c(metrics_evol$Avg_occupancy, metrics_evol$Avg_occupancy,metrics_evol$Avg_occupancy)
                         , GenModel = c(metrics_evol$GenModel,metrics_evol$GenModel,metrics_evol$GenModel)
                         , percentage_occur = c(metrics_evol$percentage_occur, metrics_evol$percentage_occur, metrics_evol$percentage_occur)
                         , n_sites = c(metrics_evol$n_sites, metrics_evol$n_sites, metrics_evol$n_sites)
                         , variable = rep(c("TSS", "Sensitivity", "Specificity"), each=length(metrics_evol$n_sites))
                         , value = c(metrics_evol$TSS_vOracle, metrics_evol$sensitivity_vOracle, metrics_evol$specificity_vOracle))
  data_all <- data_all[data_all$percentage_occur %in% c(15, 50, 80),]
  data_all <- data_all %>%
    group_by(percentage_occur, GenModel, n_sites, variable) %>%
    summarise(mean_value = mean(value, na.rm=TRUE)
              , sd_value=sd(value, na.rm=TRUE))%>%
    mutate(lower = mean_value - 2*sd_value,
           upper = mean_value + 2*sd_value)
  data_all$variable <- factor(data_all$variable, levels=c("TSS", "Sensitivity", "Specificity"))
  
  data_all$DumbModel <- gsub("Mod.", "", data_all$GenModel)
  data_all$DumbModel[data_all$DumbModel=="L"] <- "Latent"
  data_all$DumbModel <- factor(data_all$DumbModel, levels=c("Benchmark", "2V.high", "1V.high", "Latent"))
  
  evol <- ggplot(data_all, aes(x=n_sites, y=mean_value, color=factor(DumbModel), group=DumbModel))+
    geom_line()+
    geom_pointrange(aes(ymin=lower, ymax=upper, x=n_sites, y=mean_value, color=factor(DumbModel)))+
    facet_grid(variable~percentage_occur, scales = "free"
               , labeller = as_labeller(c("15"="Low occurrence"
                                          , "50"="Medium occurrence"
                                          , "80"="High occurrence"
                                          , "TSS"="Ratio TSS"
                                          , "Sensitivity"="Ratio Sensitivity"
                                          , "Specificity"="Ratio Specificity")))+
    theme_bw() +
    theme(panel.background = element_blank(),legend.position = "bottom") +
    scale_color_brewer(palette = "Set1")+
    labs(x = "Number of sites",
         y = "Average value",
         color = "Model") 
  print(evol)
  ggsave(filename = "Figures/SupFig3.jpeg", width = 7.5, height=7.5, dpi=450, units="in")
  
}

# Plots for abundance ----

# function to plot the Species Abundance Distribution (SAD)
plot_SAD <- function(data_all # formatted data from format_metrics function
                     ) {
  # data_all = metrics_data
  data_abund <- data_all[data_all$GenModel=="Benchmark",]
  # data.table::setDT(data_abund)[,c("Species_rank_PU"):=rank(Avg_Abund.tot), by=.(Universe)]
  mean_log_abund <- mean(log(data_abund$Avg_Abund.mean))
  sd_log_abund <- sd(log(data_abund$Avg_Abund.mean))
  jpeg("Figures/Fig2.jpeg", width = 6, height=4, res=300, units="in")
  
  plot <- ggplot(data_abund) +
    theme(legend.position="none") +
    geom_density(aes(x=Avg_Abund.mean, group=Universe), color="grey") +
    geom_density(aes(x=Avg_Abund.mean), stat="density")+
    geom_line(aes(x= Avg_Abund.mean, y= dlnorm(Avg_Abund.mean, mean_log_abund, sd_log_abund)),col="red")+
    xlab("Abundance") + ylab("Density") +
    theme_bw()
  gridExtra::grid.arrange(plot)
  dev.off()
  
}

# plotting the ratio and delta metrics for a fixed number of sites
plot_abund_metrics <- function(metrics_data # formatted data from format_metrics
                               ) {
  # formatting for the plots
  bins_temp <- quantile(metrics_data$Avg_occupancy, probs = seq(0, 1, 0.05), digits=7)
  metrics_data$bins_occur <- cut(metrics_data$Avg_occupancy, breaks=bins_temp, include.lowest = TRUE)
  bins_temp <- quantile(metrics_data$Avg_Abund.mean, probs = seq(0, 1, 0.05), digits=7)
  metrics_data$bins_abund <- cut(metrics_data$Avg_Abund.mean, breaks=bins_temp, include.lowest = TRUE)
  metricsP_data <- metrics_data %>%
    dplyr::group_by(GenModel, bins_occur, bins_abund) %>%
    dplyr::summarise(TSS_vOracle = mean(TSS_vOracle)
                     , sensitivity_vOracle = mean(sensitivity_vOracle)
                     , specificity_vOracle = mean(specificity_vOracle)
                     , MAPE_vOracle = mean(MAPE_vOracle)
                     , RMSPE_vOracle = mean(RMSPE_vOracle)
                     , RMSE_vOracle = mean(RMSE_vOracle)
                     , SMAPE_vOracle = mean(SMAPE_vOracle)
                     , RMRPE_vOracle = mean(RMRPE_vOracle)
                     , pred.bias_vOracle = mean(pred.bias_vOracle))
  bins_df_occur <- data.frame(percentage_occur=seq(5, 100, by=5), bins_occur=levels(metricsP_data$bins_occur))
  bins_df_occur$percentage_occur <- units::set_units(bins_df_occur$percentage_occur, "%")
  metricsP_data <- dplyr::left_join(metricsP_data, bins_df_occur, by="bins_occur")
  bins_df_abund <- data.frame(percentage_abund=seq(5, 100, by=5), bins_abund=levels(metricsP_data$bins_abund))
  bins_df_abund$percentage_abund <- units::set_units(bins_df_abund$percentage_abund, "%")
  metricsP_data <- dplyr::left_join(metricsP_data, bins_df_abund, by="bins_abund")
  library(units)
  metricsP_data$Latent <- "Environment"
  metricsP_data[grep("L", metricsP_data$GenModel),"Latent"] <- "Latent"
  metricsP_data$Latent <- factor(metricsP_data$Latent, levels = c("Environment", "Latent"))
  metricsP_data$DumbModel <- gsub("L", "", metricsP_data$GenModel)
  metricsP_data$DumbModel[metricsP_data$DumbModel=="Mod."] <- "Latent"
  metricsP_data$DumbModel <- gsub("Mod.", "", metricsP_data$DumbModel)
  metricsP_data$n.var <- factor(ifelse(grepl(2, metricsP_data$DumbModel)==TRUE, 2,
                                ifelse(grepl(1, metricsP_data$DumbModel)==TRUE,1,
                                ifelse(grepl("Benchmark", metricsP_data$DumbModel)==TRUE, 3, 0)))
                                , levels = c(3,2,1,0))
  metricsP_data$DumbModel <- factor(metricsP_data$DumbModel, levels=c("Latent"
                                                                      , "1V.low", "1V.mid", "1V.high"
                                                                      , "2V.low", "2V.mid", "2V.high", "Benchmark"))
  
  data <- data.frame(Latent = c("Environment", "Latent")
                     , Intercept = c(3.5, 3.5))
  
  metricsP_data <- metricsP_data[!metricsP_data$DumbModel %in% c("Benchmark", "Latent"),]
  # Plotting MAPE
  MAPE <- ggplot(metricsP_data, aes(percentage_abund, DumbModel)) +
    geom_tile(aes(fill=MAPE_vOracle)) +
    theme_bw() +
    theme(panel.background = element_blank()) +
    xlab("Species abundance percentiles") + ylab("Model") +
    labs(fill='Ratio MAPE') +
    scale_fill_viridis_c() +
    facet_grid(Latent~., scales = "free") +
    geom_hline(data=data,aes(yintercept = Intercept))

  
  # plotting RMSPE
  RMSPE <- ggplot(metricsP_data, aes(percentage_abund, DumbModel)) +
    geom_tile(aes(fill=RMSPE_vOracle)) +
    theme_bw() +
    theme(panel.background = element_blank()) +
    xlab("Species abundance percentiles") + ylab("Model") +
    labs(fill='Ratio RMSPE') +
    scale_fill_viridis_c() +
    facet_grid(Latent~., scales = "free") +
    geom_hline(data=data,aes(yintercept = Intercept))

  # plotting RMSE
  RMSE <- ggplot(metricsP_data, aes(percentage_abund, DumbModel)) +
    geom_tile(aes(fill=RMSE_vOracle)) +
    theme_bw() +
    theme(panel.background = element_blank()) +
    xlab("Species abundance percentiles") + ylab("Model") +
    labs(fill='Ratio RMSE') +
    scale_fill_viridis_c() +
    facet_grid(Latent~., scales = "free") +
    geom_hline(data=data,aes(yintercept = Intercept))

  # plotting SMAPE
  SMAPE <- ggplot(metricsP_data, aes(percentage_abund, DumbModel)) +
    geom_tile(aes(fill=SMAPE_vOracle)) +
    theme_bw() +
    theme(panel.background = element_blank()) +
    xlab("Species abundance percentiles") + ylab("Model") +
    labs(fill='Ratio SMAPE') +
    scale_fill_viridis_c() +
    facet_grid(Latent~., scales = "free") +
    geom_hline(data=data,aes(yintercept = Intercept))

  # plotting RMRPE
  RMRPE <- ggplot(metricsP_data, aes(percentage_abund, DumbModel)) +
    geom_tile(aes(fill=RMRPE_vOracle)) +
    theme_bw() +
    theme(panel.background = element_blank()) +
    xlab("Species abundance percentiles") + ylab("Model") +
    labs(fill='Ratio RMRPE') +
    scale_fill_viridis_c() +
    facet_grid(Latent~., scales = "free") +
    geom_hline(data=data,aes(yintercept = Intercept))
  
  
  # new column on whether the model has latence or not
  metrics_data$type <- NA
  metrics_data$type[-grep("L", metrics_data$GenModel)] <- "env"
  metrics_data$type[is.na(metrics_data$type)] <- "latent"
  # new column on a general model so we can contrast based on the previous column
  metrics_data$DumbModel <- gsub("L", "", metrics_data$GenModel)
  metrics_data$DumbModel[metrics_data$DumbModel=="Mod."] <- "Benchmark"
  metrics_data <- metrics_data[metrics_data$DumbModel!="Benchmark",]
  compar_data <- metrics_data %>% 
    dplyr::group_by(Universe, Species, DumbModel) %>%
    summarise(delta_TSS = mean(TSS[type=="env"] - TSS[type=="latent"])
              , delta_sensitivity = mean(sensitivity[type=="env"] - sensitivity[type=="latent"])
              , delta_specificity = mean(specificity[type=="env"] - specificity[type=="latent"])
              , delta_MAPE = mean(MAPE[type=="env"] - MAPE[type=="latent"])
              , delta_RMSPE = mean(RMSPE[type=="env"] - RMSPE[type=="latent"])
              , delta_RMSE = mean(RMSE[type=="env"] - RMSE[type=="latent"])
              , delta_SMAPE = mean(SMAPE[type=="env"] - SMAPE[type=="latent"])
              , delta_RMRPE = mean(RMRPE[type=="env"] - RMRPE[type=="latent"])
              , delta_pred.bias = mean(abs(pred.bias[type=="env"]) - abs(pred.bias[type=="latent"]))
              , Avg_occupancy = unique(Avg_occupancy)
              , bins_occur = unique(bins_occur)
              , bins_abund = unique(bins_abund))
  comparP_data <- compar_data %>%
    dplyr::group_by(DumbModel, bins_abund, bins_occur) %>%
    dplyr::summarise(delta_TSS = mean(delta_TSS)
                     , delta_sensitivity = mean(delta_sensitivity)
                     , delta_specificity = mean(delta_specificity)
                     , delta_MAPE = mean(delta_MAPE)
                     , delta_RMSPE = mean(delta_RMSPE)
                     , delta_RMSE = mean(delta_RMSE)
                     , delta_SMAPE = mean(delta_SMAPE)
                     , delta_RMRPE = mean(delta_RMRPE)
                     , delta_pred.bias = mean(delta_pred.bias))
  bins_df_occur <- data.frame(percentage_occur=seq(5, 100, by=5), bins_occur=levels(comparP_data$bins_occur))
  bins_df_occur$percentage_occur <- units::set_units(bins_df_occur$percentage_occur, "%")
  comparP_data <- dplyr::left_join(comparP_data, bins_df_occur, by="bins_occur")
  bins_df_abund <- data.frame(percentage_abund=seq(5, 100, by=5), bins_abund=levels(comparP_data$bins_abund))
  bins_df_abund$percentage <- units::set_units(bins_df_abund$percentage_abund, "%")
  comparP_data <- dplyr::left_join(comparP_data, bins_df_abund, by="bins_abund")
  
  comparP_data$DumbModel <- gsub("Mod.", "", comparP_data$DumbModel)
  comparP_data$DumbModel <- factor(comparP_data$DumbModel, levels=c("1V.low", "1V.mid", "1V.high"
                                                                    , "2V.low", "2V.mid", "2V.high"))
  
  MAPE_delta <- ggplot(comparP_data, aes(percentage_abund, DumbModel)) +
    geom_tile(aes(fill=delta_MAPE)) +
    theme_bw() +
    theme(panel.background = element_blank()) +
    xlab("Species abundance percentiles") + ylab("Model") +
    labs(fill="\u0394 MAPE") +
    scale_fill_gradient2(low = "blue", mid="white", high = "red", oob=scales::squish) +
    geom_hline(data=data,aes(yintercept = 3.5))

  RMSPE_delta <- ggplot(comparP_data, aes(percentage_abund, DumbModel)) +
    geom_tile(aes(fill=delta_RMSPE)) +
    theme_bw() +
    theme(panel.background = element_blank()) +
    xlab("Species abundance percentiles") + ylab("Model") +
    labs(fill='\u0394 RMSPE') +
    scale_fill_gradient2(low = "blue", mid="white", high = "red", oob=scales::squish) +
    geom_hline(data=data,aes(yintercept = 3.5))

  RMSE_delta <- ggplot(comparP_data, aes(percentage_abund, DumbModel)) +
    geom_tile(aes(fill=delta_RMSE)) +
    theme_bw() +
    theme(panel.background = element_blank()) +
    xlab("Species abundance percentiles") + ylab("Model") +
    labs(fill='\u0394 RMSE') +
    scale_fill_gradient2(low = "blue", mid="white", high = "red", oob=scales::squish) +
    geom_hline(data=data,aes(yintercept = 3.5))
  
  SMAPE_delta <- ggplot(comparP_data, aes(percentage_abund, DumbModel)) +
    geom_tile(aes(fill=delta_SMAPE)) +
    theme_bw() +
    theme(panel.background = element_blank()) +
    xlab("Species abundance percentiles") + ylab("Model") +
    labs(fill='\u0394 SMAPE') +
    scale_fill_gradient2(low = "blue", mid="white", high = "red", oob=scales::squish) +
    geom_hline(data=data,aes(yintercept = 3.5))

  RMRPE_delta <- ggplot(comparP_data, aes(percentage_abund, DumbModel)) +
    geom_tile(aes(fill=delta_RMRPE)) +
    theme_bw() +
    theme(panel.background = element_blank()) +
    xlab("Species abundance percentiles") + ylab("Model") +
    labs(fill='\u0394 RMRPE') +
    scale_fill_gradient2(low = "blue", mid="white", high = "red", oob=scales::squish) +
    geom_hline(data=data,aes(yintercept = 3.5))

  jpeg("Figures/SupFig4.jpeg", width = 8, height=12, res=600, units="in")
  
  gridExtra::grid.arrange(MAPE, RMSPE, RMSE, SMAPE, RMRPE
                          , MAPE_delta, RMSPE_delta, RMSE_delta, SMAPE_delta, RMRPE_delta
                          , ncol=2, as.table=FALSE)
  dev.off()
  
  cowplot::plot_grid(MAPE, MAPE_delta, nrow=2, align = "v", axis = "Model", rel_heights = c(1.5,1))
  # print(MAPE)
  ggsave(filename = "Figures/Fig6.jpeg", width = 8, height=8, dpi=300, units="in")
  
}

# plotting the correlation of metrics depending on models studied
plot_abund_correlation <- function(metrics_data # formatted data from format_metrics
                                   ) {
  bins_temp <- quantile(metrics_data$Avg_Abund.mean, probs = seq(0, 1, 0.05), digits=7)
  metrics_data$bins_abund <- cut(metrics_data$Avg_Abund.mean, breaks=bins_temp, include.lowest = TRUE)
  bins_temp <- quantile(metrics_data$Avg_occupancy, probs = seq(0, 1, 0.05), digits=7)
  metrics_data$bins_occur <- cut(metrics_data$Avg_occupancy, breaks=bins_temp, include.lowest = TRUE)
  
  data_oracle <- format_oracle(n.species, universes, n.sites)
  bins_temp <- quantile(data_oracle$Avg_Abund.mean, probs = seq(0, 1, 0.05), digits=7)
  data_oracle$bins_abund<- cut(data_oracle$Avg_Abund.mean, breaks=bins_temp, include.lowest = TRUE)
  bins_temp <- quantile(data_oracle$Avg_occupancy, probs = seq(0, 1, 0.05), digits=7)
  data_oracle$bins_occur <- cut(data_oracle$Avg_occupancy, breaks=bins_temp, include.lowest = TRUE)
  
  metrics_data <- metrics_data[metrics_data$GenModel %in% c("Benchmark", "Mod.L"), ]
  metrics_data <- metrics_data[,colnames(metrics_data) %in% colnames(data_oracle)]
  all_data <- rbind(metrics_data, data_oracle)
  bins_df_abund <- data.frame(percentage_abund=seq(5, 100, by=5), bins_abund=levels(all_data$bins_abund))
  bins_df_occur <- data.frame(percentage_occur=seq(5, 100, by=5), bins_occur=levels(all_data$bins_occur))
  all_data <- dplyr::left_join(all_data, bins_df_abund, by="bins_abund")
  all_data <- dplyr::left_join(all_data, bins_df_occur, by="bins_occur")
  all_data$Model <- factor(ifelse(all_data$GenModel=="Mod.L", "Latent", all_data$GenModel)
                           , levels = c("Oracle", "Benchmark", "Latent"))
  data_PA <- data.frame(Variable = rep(c("MAPE", "RMSPE", "RMSE", "SMAPE", "RMRPE"), each = length(all_data$RMRPE))
                        , Value = c(all_data$MAPE, all_data$RMSPE, all_data$RMSE, all_data$SMAPE, all_data$RMRPE)
                        , percentage_abund = c(all_data$percentage_abund, all_data$percentage_abund, all_data$percentage_abund, all_data$percentage_abund, all_data$percentage_abund)
                        , Model = c(all_data$Model, all_data$Model, all_data$Model, all_data$Model, all_data$Model))
  data_PA <- data_PA %>%
    group_by(Variable, percentage_abund, Model) %>%
    summarise(mean_value = mean(Value),
              sd_value = sd(Value,na.rm= TRUE))%>%
    mutate(lower = mean_value - 2*sd_value,
           upper = mean_value + 2*sd_value)
  data_PA$Variable <- factor(data_PA$Variable, levels=c("MAPE", "RMSPE", "RMSE", "SMAPE", "RMRPE"))
  jpeg("Figures/SupFig5.jpeg", width = 5, height=6, res=300, units="in")
  
  TSS <- ggplot(data_PA, aes(x=percentage_abund, y=mean_value, group=Model, color=Model))+
    geom_line() +
    geom_pointrange(data = data_PA, aes(ymin=lower, ymax=upper, x=percentage_abund, y=mean_value, color=Model))+
    facet_grid(Variable~., scales = "free_y")+
    theme_bw() +
    theme(panel.background = element_blank(),legend.position = "bottom") +
    scale_color_brewer(palette = "Set1")+
    xlab("Species abundance percentiles") + ylab("Average value") 
  gridExtra::grid.arrange(TSS)
  dev.off()
}

# plotting the average metric depending on number of sites used to fit the data, number of species
plot_abund_sites <- function(metrics_evol # formatted data from format_evol
                             ) {
  metrics_evol <- metrics_evol[metrics_evol$GenModel %in% c("Benchmark", "Mod.L", "Mod.1V.high", "Mod.2V.high"),]
  metrics_evol <- metrics_evol[,colnames(metrics_evol) %in% c("Species", "MAPE_vOracle", "RMSPE_vOracle"
                                                              , "RMSE_vOracle", "SMAPE_vOracle", 
                                                              "RMRPE_vOracle", "GenModel"
                                                              , "Avg_Abund.tot", "Universe", "n_sites")]
  bins_temp <- quantile(metrics_evol$Avg_Abund.tot, probs = seq(0, 1, 0.05), digits=7)
  metrics_evol$bins_abund <- cut(metrics_evol$Avg_Abund.tot, breaks=bins_temp, include.lowest = TRUE)
  bins_df_abund <- data.frame(percentage_abund=seq(5, 100, by=5), bins_abund=levels(metrics_evol$bins_abund))
  metrics_evol <- dplyr::left_join(metrics_evol, bins_df_abund, by="bins_abund")
  
  data_all <- data.frame(Universe = c(metrics_evol$Universe, metrics_evol$Universe, metrics_evol$Universe,  metrics_evol$Universe,  metrics_evol$Universe)
                         , Avg_Abund.tot = c(metrics_evol$Avg_Abund.tot, metrics_evol$Avg_Abund.tot,metrics_evol$Avg_Abund.tot,metrics_evol$Avg_Abund.tot,metrics_evol$Avg_Abund.tot)
                         , GenModel = c(metrics_evol$GenModel,metrics_evol$GenModel,metrics_evol$GenModel,metrics_evol$GenModel,metrics_evol$GenModel)
                         , percentage_abund = c(metrics_evol$percentage_abund, metrics_evol$percentage_abund, metrics_evol$percentage_abund, metrics_evol$percentage_abund, metrics_evol$percentage_abund)
                         , n_sites = c(metrics_evol$n_sites, metrics_evol$n_sites, metrics_evol$n_sites, metrics_evol$n_sites, metrics_evol$n_sites)
                         , variable = rep(c("MAPE", "RMSPE", "RMSE", "SMAPE", "RMRPE"), each=length(metrics_evol$n_sites))
                         , value = c(metrics_evol$MAPE_vOracle, metrics_evol$RMSPE, metrics_evol$RMSE, metrics_evol$SMAPE_vOracle, metrics_evol$RMRPE_vOracle))
  data_all <- data_all[data_all$percentage_abund %in% c(15, 50, 80),]
  data_all <- data_all %>%
    group_by(percentage_abund, GenModel, n_sites, variable) %>%
    summarise(mean_value = mean(value, na.rm=TRUE)
              , sd_value=sd(value, na.rm=TRUE))%>%
    mutate(lower = mean_value - 2*sd_value,
           upper = mean_value + 2*sd_value)
  data_all$variable <- factor(data_all$variable, levels=c("MAPE", "RMSPE", "RMSE", "SMAPE", "RMRPE"))
  
  
  data_all$DumbModel <- gsub("Mod.", "", data_all$GenModel)
  data_all$DumbModel[data_all$DumbModel=="L"] <- "Latent"
  data_all$DumbModel <- factor(data_all$DumbModel, levels=c("Benchmark", "2V.high", "1V.high", "Latent"))
  
  jpeg("Figures/SupFig6.jpeg", width = 7.5, height=9, res=300, units="in")
  
  evol <- ggplot(data_all, aes(x=n_sites, y=mean_value, color=factor(DumbModel), group=DumbModel))+
    geom_line()+
    geom_pointrange(aes(ymin=lower, ymax=upper, x=n_sites, y=mean_value, color=factor(DumbModel)))+
    facet_grid(variable~percentage_abund, scales = "free"
               , labeller = as_labeller(c("15"="Low abundance"
                                          , "50"="Medium abundance"
                                          , "80"="High abundance"
                                          , "MAPE"="Ratio MAPE"
                                          , "RMSPE"="Ratio RMSPE"
                                          , "RMSE"="Ratio RMSE"
                                          , "SMAPE"="Ratio SMAPE"
                                          , "RMRPE"="Ratio RMRPE")))+
    theme_bw() +
    theme(panel.background = element_blank(),legend.position = "bottom") +
    scale_color_brewer(palette = "Set1")+
    labs(x = "Number of sites",
         y = "Average value",
         color = "Model") 
  print(evol)
  dev.off()
}
