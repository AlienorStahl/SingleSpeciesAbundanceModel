# Evaluating number of latents ----

# Init ----
  
  # Libraries ----
  library(ecoCopula) #needed for latent var analysis
  library(vegan)     #needed for RDA analysis
  library(tidyr)     #used for setting up sim data
  library(ggformula) #used for plotting
  library(dplyr)     #needed for data munging
  library(progress)  #create a progress bar
  

  # Functions ----
  #function to create new regression coefficients for each landscape
  generate_covs <-  function(n_sp # number of species
                             , n_dim # number of dimensions (environmental variables)
                             , intercept_range = c(-2.4,1.2) # range of values for the intercepts
                             , cov_range = c(-0.8, 0.8) # range of values for the slopes
                             ){
    intercept <- runif(n_sp, intercept_range[1], intercept_range[2])
    slopes <- matrix(runif(n_sp*n_dim,-0.8,0.8),nrow = n_dim,ncol = n_sp)
    bX <- rbind(intercept, slopes )
    colnames(bX) <- paste0("Spc", 1:n_sp)
    rownames(bX) <- paste0("b", 0:n_dim)
    bX
  }
  
  #function to generate grid of abundances
  generate_community <- function(n_sites # number of sites in the landscape, 
                                 , n_sp # number of species in the community
                                 , n_dim # number of independent environmental variables
                                 , bX # matrix of the environmental variables
                                 ) {
    X <- cbind(rep(1, times=n_sites),
               matrix(data = rnorm(n_sites * n_dim, mean=0, sd=1),
                      nrow = n_sites,
                      ncol = n_dim))
    colnames(X) <- c("Int", paste0("X", 1:n_dim))
    Y <- exp(X %*% bX)
    Y <- rpois(n_sites*n_sp,Y)
    Y <- matrix(Y,n_sites,n_sp,byrow=FALSE)
    colnames(Y) <- paste0("Sp", 1:n_sp) 
    return(list(Y=Y, X=X))
  }
  
  
  
  
# I. Simulations ----
  #create a grid of predictor variables.
  #varying # of species, dimensions of the predictor, number of sites in the universe
  # and have 10 replicates for each combo
  input_data <- expand_grid(n_sp =  c(10, 20, 30) 
                            ,n_dim = 1:5 
                            ,n_sites = c(100,200, 300)
                            ,rep = 1:10) %>%
    mutate(seed = 1:n())
  
  
  n_sims <- nrow(input_data)
  
  #specify how many latent variables to try
  max_latents <- 5
  
  output_data <- expand_grid(input_data, n_latents = 1:max_latents)%>%
    mutate(var_exp = NA,
           BIC = NA)
  
  #keeps track of where to send output
  counter <- 1
  #creates a progress bar
  pb <- progress_bar$new(total = n_sims)
  
  pb$tick(0)
  for(i in 1:n_sims){
    # setting the seed
    set.seed(input_data$seed[i])
    # parameters of the current simulation
    n_dim <- input_data$n_dim[i]
    n_sp <- input_data$n_sp[i]
    n_sites <- input_data$n_sites[i]
    # generating the environmnetal variables
    covs <- generate_covs(n_sp, n_dim)
    # generating the communities
    landscape <- generate_community(n_sites, n_sp, n_dim, covs)
    
    # formatting abundance, environment and PA
    counts <- landscape$Y
    X <- landscape$X[,-1] #remove the intercept
    pa <- (counts>0) + 0
    
    #fits a null (intercept-only) SDM model for each species
    pa_sdm <- stackedsdm(pa, ~1, data = X, family = "binomial")
    
    #look at each number of latents in turn.
    for(j in 1:max_latents){
      
      #the try() here is because rarely factanal returns an error for
      #the randomly generated starting coordinates; this tries a second
      #time, then returns NA if that doesn't work
      latents <- try(cord(pa_sdm, j),silent = TRUE)
      if(class(latents)=="try-error"){
        #retry it one last time.
        latents <- try(cord(pa_sdm, j),silent = TRUE)
      }
      if(class(latents)=='try-error'){
        output_data$var_exp[counter] <- NA
        output_data$BCI[counter] <- NA
      }else{
        #R-squared calculated w/ an RDA of latent variables predicting X
        latent_rda <- rda(scale(X) ~ scale(latents$scores))
        output_data$var_exp[counter] <- RsquareAdj(latent_rda)$adj.r.squared
        output_data$BIC[counter] <- latents$BIC
      }
      
      counter <- counter + 1
    }
    pb$tick()
  }
  # saving the results to avoid having to re-run the simulations
  saveRDS(file="Results/evaluating-latents.RData", output_data)
  # and importing them
  output_data2 <- readRDS(file="Results/evaluating-latents.RData")
  
  #calculate mean and sd of outcomes. 
  summary_data <- output_data %>%
    group_by(n_sp, n_dim, n_sites, rep)%>%
    mutate(delta_BIC = BIC - min(BIC))%>%
    group_by(n_sp, n_dim, n_latents,n_sites)%>%
    summarize(var_exp_mean = mean(var_exp,na.rm = TRUE),
              var_exp_sd = sd(var_exp,na.rm= TRUE),
              delta_BIC_mean  = mean(delta_BIC))%>%
    mutate(lower = var_exp_mean - 2*var_exp_sd,
           upper = var_exp_mean + 2*var_exp_sd)
  
# II. Plots ----
  #plotting r-squared
  plt1 <- summary_data %>%
    mutate(n_dim = paste(n_dim, "dimensional\nenvironment"),
           n_sites = paste(n_sites, "sites"))%>%
    gf_line(var_exp_mean ~ n_latents, data = ., color = ~factor(n_sp)) %>%
    gf_pointrange(var_exp_mean+ lower+upper ~ n_latents)%>%
    gf_labs(x = "Number of latent variables",
            y = "Explained inertia",
            color = "Number of species")%>%
    gf_facet_grid(n_dim~n_sites) %>%
    gf_refine(scale_color_brewer(palette = "Set1"))%>%
    gf_theme(theme_bw()) %>%
    gf_theme(legend.position = "bottom")
  
  #plotting delta-BIC
  plt2 <- summary_data %>%
    mutate(n_dim = paste(n_dim, "dimensional\nenvironment"),
           n_sites = paste(n_sites, "sites"))%>%
    gf_line(delta_BIC_mean ~ n_latents, data = ., color = ~factor(n_sp)) %>%
    gf_labs(x = "Number of latent variables",
            y = "average \u0394 BIC",
            color = "Number of species")%>%
    gf_facet_grid(n_dim~n_sites) %>%
    gf_refine(scale_color_brewer(palette = "Set1"))%>%
    gf_theme(theme_bw()) %>%
    gf_theme(legend.position = "bottom")
  
  # printing and saving the figures into a separate folder
  # all existing scenarios
  print(plt1)
  ggsave(filename = "Figures/SupFig1.jpeg", width = 8, height=9, dpi=600, units="in")
  print(plt2)
  ggsave(filename = "Figures/SupFig2.jpeg", width = 8, height=9, dpi=600, units="in")
  
  # scenario with 200 sites as a figure in the main text
  plt3 <- summary_data[summary_data$n_sites==200,] %>%
    mutate(n_dim = paste(n_dim, "dimensional\nenvironment"),
           n_sites = paste(n_sites, "sites"))%>%
    gf_line(var_exp_mean ~ n_latents, data = ., color = ~factor(n_sp)) %>%
    gf_pointrange(var_exp_mean+ lower+upper ~ n_latents)%>%
    gf_labs(x = "Number of latent variables",
            y = "Explained inertia",
            color = "Number of species")%>%
    gf_facet_grid(n_dim ~.) %>%
    gf_refine(scale_color_brewer(palette = "Set1"))%>%
    gf_theme(theme_bw()) %>%
    gf_theme(legend.position = "bottom")
  print(plt3)
  ggsave(filename = "Figures/Fig3.jpeg", width = 7, height=8, dpi=600, units="in")
  
