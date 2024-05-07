# Figures ----

# Init ----

  # Libraries ----
  library(ggplot2)    # for plots
  library(patchwork)  # to format graphs together
  library(dplyr)      # for pipelines
  library(reshape2)   # for melt

  # Functions ----
  source("Scripts/Figures_functions.R")
  
# end Init

# Figure 1 ----
  # setting the seed
  set.seed(416)
  
  # parameters of the community
  sample.size = 122
  n.species=10
  
  # environmental variables for the community
  bX <- rbind(runif(n.species, -2.4, 1.2), matrix(runif(n.species*2,-0.8,0.8),2,n.species))
  colnames(bX) <- paste0("Spc", 1:n.species)
  rownames(bX) <- c("b0", "b1", "b2")
  X <- cbind(rep(1, times=sample.size), rep(seq(-3,3, 0.1), each=2), rep(seq(-3,3,0.1), times=2))
  colnames(X) <- c("b0", "X1", "X2")
  
  # plotting the first env variable
  X1 <- ggplot() +
    geom_point(data=data.frame(X), aes(y=X1, x=1:122))+
    labs(y="X1", x="")+ 
    theme(axis.text.x =element_blank(), 
          axis.ticks.x=element_blank()) +
    theme_bw(base_size=17)
  X2 <- ggplot() +
    geom_point(data=data.frame(X), aes(y=X2, x=1:122))+
    labs(y="X2", x="Lakes")+ 
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank()) +
    theme_bw(base_size=17)
  
  
  # generating the abundance and PA
  abundance <- generate_community(sample.size, X, bX, n.species)$Y
  PA <- ifelse(abundance==0, 0, 1)
  
  # formatting it
  data_abund <- data.frame(Species = rep(paste0("Sp", 1:10), each=122)
                           , abund = (as.vector(abundance))
                           , Lakes = rep(1:122, n.species))
  data_abund$Species <- factor(data_abund$Species, level=paste0("Sp", 1:10))
  data_abund$abund[data_abund$abund==0] <- NA
  
  # plotting the abundance of each species
  abund_OC <- ggplot() +
    geom_tile(data=data_abund, aes(y=Species, x=Lakes, fill=abund)) +
    scale_fill_viridis_c(na.value="White")+  labs(x="Lakes", fill="Abundance")+ 
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) +
    theme_bw(base_size=17)
  
  
  # formatting for the PA plot 
  data_new_plot <- data.frame(Species = rep(paste0("Sp", 1:10), each=122)
                              , PA = factor(as.vector(PA))
                              , Abundance = as.vector(abundance)
                              , X1 = c(rep(X[,c("X1")], n.species), rep(X[,c("X1")], n.species), rep(X[,c("X1")], n.species))
                              , X2 = c(rep(X[,c("X2")], n.species), rep(X[,c("X2")], n.species), rep(X[,c("X2")], n.species)))
  data_new_plot$lake <- rep(1:122, n.species) 
  data_new_plot$Species <- factor(data_new_plot$Species, level=paste0("Sp", 1:10))
  
  
  # data to fit the model
  eco_lvs <- copula_lat(PA)
  data_new_plot$lat <- rep(eco_lvs, n.species)
  # plotting the latent variable compared to the environment
  Lat.plot <- ggplot()+
    geom_point(data=data_new_plot, aes(y=lat, x=lake), size=.8)+
    labs(y="Latent", x="")+ 
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank()) +
    theme_bw(base_size=17)
  
  # dataframe to record the results
  PA.pred <- abund.pred <- data.frame(matrix(data=NA, nrow=122, ncol=n.species))
  # calculating the expected PA and abundance for each species
  for (i in 1:n.species) {
    species.temp <- paste0("Sp", i)
    data_all <- data.frame(cbind(abundance[,species.temp], eco_lvs, X[,c("X1")]))
    colnames(data_all) <- c("Sp", "Factor1", "X1")
    mod_lat <- mvabund::manyglm(formula = Sp ~ X1 + Factor1, family="poisson", data=data_all) 
    predictions_lat <- predict(mod_lat, type="response")
    # PA predictions
    pred.presence.lat = 1-ppois(0, predictions_lat)
    threshold.species <-sum(PA[,i]/sample.size)
    lat.PA <- ifelse(pred.presence.lat > threshold.species, 1, 0)
    PA.pred[,i] <- lat.PA
    # Abundance
    abundance.lat <- predictions_lat/pred.presence.lat
    abundance.lat[PA[,i]==0] <- NA
    abund.pred[,i] <- abundance.lat
  }
  
  # formatting the data for future plots
  data_pred <- data.frame(Species = rep(paste0("Sp", 1:10), each=122)
                          , PA = factor(unlist(PA.pred))
                          , Abundance = unlist(abund.pred)
                          , Lakes = rep(1:122, n.species))
  data_pred$Species <- factor(data_pred$Species, level=paste0("Sp", 1:10))
  
  # predicting presence absence for a single species
  PA_lat <- ggplot() +
    geom_tile(data=data_pred[data_pred$Species=="Sp10",], aes(x=Lakes, y=Species, fill=PA)) +
    scale_fill_manual(labels = c("0"="Absent", "1"="Present"), values=c("white", "black"))+
    labs(fill="", y=expression(paste("Target \n species")), x="", title="Presence-absence")+ 
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    theme_bw(base_size=17)
  
  # predicting abundance for a single species
  Abund_lat <- ggplot() +
    geom_tile(data=data_pred[data_pred$Species=="Sp10",], aes(y=Species, x=Lakes, fill=Abundance))+
    scale_fill_viridis_c(na.value="White", limits =c(min(data_abund$abund, na.rm = TRUE), max(data_abund$abund, na.rm = TRUE)))+
    labs(fill="Abundance", y=expression(paste("Target \n species")), x="Lakes", title="Predicted abundance")+ 
    theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    theme_bw(base_size=17)
  
  # formatting for the plot
  design <- "
  15#
  15#
  156
  257
  35#
  45#
  "
  # plotting everything together and saving it
  abund_OC + Lat.plot + X1 + X2 + plot_spacer() + PA_lat + Abund_lat + plot_layout(design = design, guides = "collect", widths = c(1, 0.2, 0.8))
  ggsave(filename = "Figures/Fig1.jpeg", width = 14, height=9, dpi=600, units="in")
  
# end Fig 1

# Figures from the simulations ----
  # Parameters of interest
  n.species <- 20
  universes <- 1:30
  n.sites <- 200
  
  # importing and formating the data for a specific number of sites
  metrics_data <- format_metrics(n.species, universes, n.sites)
  # importing and formatting data for all possible scenarios
  metrics_evol <- format_evol(n.species, universes)
  
  # plotting figure 2 (species abundance distribution)
  plot_SAD(metrics_data)

  # plotting figure 4 (ratio TSS and delta TSS)
  plot_PA_metrics(metrics_data)
  # plotting figure 5 (correlation of TSS)
  plot_PA_correlation(metrics_data)
  # plotting sup figure 3
  plot_PA_sites(metrics_evol)
  
  # plotting fig 6 and sup fig 4
  plot_abund_metrics(metrics_data)
  # plotting sup fig 5
  plot_abund_correlation(metrics_data)
  # plotting sup fig 6
  plot_abund_sites(metrics_evol)
  
  
# end ----
  