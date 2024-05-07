# Project on Advancing single species abundance models: robust models for predicting abundance using co-occurrence from communities

Author names: Aliénor Stahl (1), Eric J. Pedersen (1) and Pedro R. Peres-Neto (1)
Affiliations: 1: Concordia University, Department of Biology, Montreal, Canada

Scripts are within the "Scripts" folder. All figures generated are saved in the "Figures" folder. Results generated throughout the simulations are saved in the "Results" folder.

# Scripts
evaluating-latents: Script used to evaluate the number of latent variables needed to capture environmental variation. It starts by simulating communites and calculating the BIC and R squared as functions of environmental variables, number of species, number of sites, and number of latent variables. The results are saved in the "evaluating-latents.RData" before being plotted
Simulations: Script used to assess the performance of species abundance models. Models differ in the predictors included. Several scenarios of sample size used to fit the models are used. The results are saved per landscape and within files referring to the sample size (i.e., Sample.size_100). It uses functions defined in the Simulations_functions script. The results for each landscape are saved within a Sample.size_X folder, with X referring to the number of sites used to fit the various models within each landscape. Each landscape is saved as a different Rdata file.
Simulations_functions: Functions used in the Simulations script. Functions range from generating communities, fitting models, generating latent variables to generating metrics of errors based on the models' predictions
Figures: Script used to generate the final figures of the article. It uses functions defined in the Figures_functions script.
Figures_functions: Functions used in the Figures script.

# Results
evaluating-latents.RData: Generated from the evaluating-latents script
Sample.size_X folders: X refers to the number of sites used to fit the various models within each landscape (X = [100, 200, 300, 500]). Each folder contains 30 RData files named as "Simul_universe_Y" with Y refering to the landscape. (Y takes values of integer between 1 and 30). These data are generated from the Simulations script and are used to plot the figures in the Figures script.

# Figures
Figures generated from the Simulations and Figures scripts. 
