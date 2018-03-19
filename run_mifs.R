#######################################################################################################
# Script example of maximum iterated filtering
# The search was initiated from 1e2 starting parameter sets, sampled within the confidence intervals
# of the trajectory matching
# To save time, this script runs a single estimation, with 2 MIF iterations
# CAUTION: the data used in this script are NOT the real data, but synthetic data simulated with the best model
#######################################################################################################

rm(list=ls())

# Loads the required packages, checks the version number ------------------
library(checkpoint)
checkpoint(snapshotDate = "2015-10-01", scanForPackages = F)

library(nloptr)
library(magrittr)
library(plyr)
library(reshape2)
library(ggplot2)
library(subplex)
library(pomp)
print(packageVersion("pomp")) # Should be 1.2.1.1
source("create_pomp_model.R") # Source function that creates the pomp object
model <- "WV" # Model to be estimated; SIR: no loss; WV: waning; LV: leaky vaccine; LV_WV: leaky+waning model
theme_set(theme_bw())

# Create age categories and data ------------------------------------------
# Age classes in the model
agecats_mod <- c("[0,0.33)", "[0.33,1)", "[1,5)", "[5,10)", "[10,15)",
                 "[15,20)", "[20,25)", "[25,30)", "[30,35)", "[35,40)", "[40,45)",
                 "[45,50)", "[50,55)", "[55,60)" , "[60,65)" , "[65,70)" , "[70,75)")
nages_mod <- length(agecats_mod) # 17 age classes

# Load monthly case reports, 1990--2005 (16 * 12 = 192 time points)
# These are simulated data from the waning model
data_month <- readRDS("simulated_data.rds")

# Fix model parameters --------------------------------------------------------
nbasis <- 3L # No of periodic spline bases for seasonal functions
time_start_sim <- -120 # Time at which the simulations should start (relative to 1990)
time_start_data <- 1990 # Time at which the data start
covars <- readRDS("covars.rds") # Covariates

# Append 'fake' NA observation before the first data point
time_append <- seq(time_start_sim + 1, 0, by = 1) + 1990
data_append <- matrix(data = NA, nrow = length(time_append), 
                      ncol = ncol(data_month), 
                      dimnames = list(NULL, colnames(data_month))) %>%
  as.data.frame(.) %>% mutate(time = time_append)
data_month <- rbind(data_append, data_month)

# Create pomp object
# Time step of 0.01 yr
VSEIR2 <- create_pomp_model(data = data_month, 
                            covars = covars, 
                            nages_mod = nages_mod, 
                            nbasis = nbasis,
                            time_start_data = time_start_data, 
                            time_start_sim = time_start_sim, 
                            dt = 0.01)
theta <- readRDS("mle_waning_model_stochastic.rds") # MLE of the waning model
coef(VSEIR2) <- theta

# List of estimated parameters, depending on the model --------------------
switch(model,
       # No-loss model
       SIR = {
         coef(VSEIR2, c("alphaV", "epsilonV", "theta", "rho2")) <- 0
         estpars <- c("q1-1", "q1-2", "q1-3", 
                      paste0("omegaC", 1:(nbasis - 1)), 
                      paste0("omegaT", 1:(nbasis - 1)), 
                      "rho1-4", "tau", "epsA")
       },
       # Waning model: waning vaccinal immunity
       WV = {
         coef(VSEIR2, c("epsilonV")) <- 0
         estpars <- c("q1-1", "q1-2", "q1-3", 
                      paste0("omegaC", 1:(nbasis - 1)), 
                      paste0("omegaT", 1:(nbasis - 1)), 
                      "alphaV", "rho1-4", "tau", "rho2", "theta", "epsA")
       },
       # Leaky model: leaky vaccinal immunity
       LV = {
         coef(VSEIR2, c("alphaV")) <- 0
         estpars <- c("q1-1", "q1-2", "q1-3", 
                      paste0("omegaC", 1:(nbasis - 1)), 
                      paste0("omegaT", 1:(nbasis - 1)), 
                      "epsilonV", "rho1-4", "tau", "rho2", "theta", "epsA")
       },
       # Leaky + waning model: waning and leaky vaccinal immunity
       LV_WV = {
         estpars <- c("q1-1", "q1-2", "q1-3", 
                      paste0("omegaC", 1:(nbasis - 1)), 
                      paste0("omegaT", 1:(nbasis - 1)), 
                      "epsilonV", "alphaV", "rho1-4", "tau", "rho2", "theta", "epsA")
       }
)

# Start range for estimated parameters ------------------------------------
q1_range <- matrix(c(0.05, 0.1, 0.3, 0.7, 0.1, 0.2), ncol = 3, nrow = 2) %>% as.data.frame(.)
omegaC_range <- matrix(c(-1.4, 0.5, 0.7, 1.3), ncol = 2, nrow = 2) %>% as.data.frame(.)
omegaT_range <- matrix(c(-0.5, 0.2, -0.8, 0.1), ncol = 2, nrow = 2) %>% as.data.frame(.)
rho1_range <- matrix(c(0.23, 1), ncol = 1, nrow = 2) %>% as.data.frame(.)
start.range <- data.frame(epsA = c(0.02, 0.07),
                          q1_range, 
                          omegaC_range, 
                          omegaT_range, 
                          theta = c(0.46, 1),
                          epsilonV = c(0.03, 0.17), 
                          alphaV = c(0.007, 0.044), 
                          rho1_range, 
                          rho2 = c(0.15, 1), 
                          tau = c(0.50, 0.62))
colnames(start.range) <- c("epsA", 
                           paste0("q1-", 1:3), 
                           paste0("omegaC", 1:2), 
                           paste0("omegaT", 1:2),
                           "theta", 
                           "epsilonV", 
                           "alphaV",
                           "rho1-4", 
                           "rho2", 
                           "tau")
start.range <- subset(start.range, select = estpars)
start.values <- sobolDesign(lower = setNames(as.numeric(start.range[1, ]), names(start.range[1, ])), 
                            upper = setNames(as.numeric(start.range[2, ]), names(start.range[2, ])), 
                            nseq = 1e2)

# Run MIFs -------------------------------------------------

# Set the intensity of the random walks for all the parameters
ran_sd1 <- ran_sd2 <- list()
ran_sd_val <- 1e-2 # Intensity of random walk
for(i in seq_along(estpars)) {
  ran_sd1[[estpars[i]]] <- expression(ifelse(time > 0, 1e-6, 0))
  ran_sd2[[estpars[i]]] <- expression(ifelse(time > 0, ran_sd_val, 0))
}

for(i in 1:1) {
  cat("Estimation no ", i, "\n")
  x0 <- as.numeric(start.values[i, estpars]) # Starting parameters on the natural scale
  coef(VSEIR2, estpars) <- x0
  
  # Run MIF
  # This example runs 2 MIF iterations with 2e2 particles
  # The real results were obtained with 50 MIF iterations and 2e3 particles
  seed <- 21862186L
  tic <- Sys.time()
  mf <- freeze(seed = seed, expr = {
    mif2(object = VSEIR2, 
         Nmif = 1, 
         Np = 2e2, 
         start = coef(VSEIR2),
         cooling.type = "hyperbolic", 
         cooling.fraction.50 = 0.1,
         rw.sd = ran_sd1, 
         transform = TRUE, 
         verbose = T) %>% 
      continue(Nmif = 1, rw.sd = ran_sd2, verbose = T)
  })
  toc <- Sys.time()
  print(toc - tic)
  
  # Calculate the likelihood at the estimated parameter set
  # This example runs 2 replicate particle filters with 2e2 particles
  # The results were based on 20 replicate particle filters 5e4 particles
  coef(VSEIR2, estpars) <- unname(coef(mf, estpars)) # Set the values to the parameter estimates
  tic <- Sys.time()
  llpost <- freeze(seed = seed, expr = {
    replicate(2, logLik(pfilter(VSEIR2,
                                Np = 2e2,
                                verbose = T,
                                params = coef(VSEIR2),
                                max.fail = Inf)))
  })
  toc <- Sys.time()
  print(toc - tic)
  
  # Print mean log-likelihood and associated SE
  print(logmeanexp(llpost, se = T))
}

#######################################################################################################
# End
#######################################################################################################