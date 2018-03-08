#######################################################################################################
# Script example of trajectory matching 
# The broad search was initiated from 1e4 starting parameter sets, with maximum execution time of 20 min
# for each estimation
# To save time, this script runs a single estimation, with maximum execution time of 1 min
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
model <- "WV" # Model to be estimated; SIR: no loss; WV: waning; LV: leaky vaccine

# Create age categories and data ------------------------------------------
# Age classes in the model
agecats_mod <- c("[0,0.33)", "[0.33,1)", "[1,5)", "[5,10)", "[10,15)",
                 "[15,20)", "[20,25)", "[25,30)", "[30,35)", "[35,40)", "[40,45)",
                 "[45,50)", "[50,55)", "[55,60)" , "[60,65)" , "[65,70)" , "[70,75)")
nages_mod <- length(agecats_mod) # 17 age classes

# Load data
data_month <- readRDS("data.rds")

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
VSEIR2 <- create_pomp_model(data = data_month, 
                            covars = covars, 
                            nages_mod = nages_mod, 
                            nbasis = nbasis,
                            time_start_data = time_start_data, 
                            time_start_sim = time_start_sim, 
                            dt = 0.01)
theta <- readRDS("mle_waning_model_deterministic.rds")
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
q1_range <- matrix(c(0, 1), ncol = 3, nrow = 2) %>% as.data.frame(.)
omegaC_range <- matrix(c(-5, 5), ncol = (nbasis - 1), nrow = 2) %>% as.data.frame(.)
omegaT_range <- matrix(c(-5, 5), ncol = (nbasis - 1), nrow = 2) %>% as.data.frame(.)
rho1_range <- matrix(c(0, 1), ncol = 1, nrow = 2) %>% as.data.frame(.)
start.range <- data.frame(epsA = c(0, 1),
                          q1_range, omegaC_range, omegaT_range, 
                          theta = c(0, 1),
                          epsilonV = c(0, 1), 
                          alphaV = c(0, 10), 
                          rho1_range, 
                          rho2 = c(0, 1), 
                          tau = c(0, 10))
colnames(start.range) <- c("epsA", 
                           paste0("q1-", 1:3), 
                           paste0("omegaC", 1:(nbasis - 1)), 
                           paste0("omegaT", 1:(nbasis - 1)),
                           "theta", 
                           "epsilonV", 
                           "alphaV",
                           "rho1-4", 
                           "rho2", 
                           "tau")
start.range <- subset(start.range, select = estpars)
start.values <- sobolDesign(lower = setNames(as.numeric(start.range[1, ]), names(start.range[1, ])), 
                            upper = setNames(as.numeric(start.range[2, ]), names(start.range[2, ])), 
                            nseq = 1e4)

# Run trajectory matching -------------------------------------------------
ofun <- traj.match.objfun(object = VSEIR2, est = estpars, transform = TRUE) # Create the objective function
for(i in 1:1){
  cat("Estimation: ", i, "\n")
  x0 <- as.numeric(start.values[i, ])
  coef(VSEIR2, estpars) <- x0
  x0_trans <- coef(VSEIR2, estpars, transform = TRUE)
  
  # Runs trajectory matching using the subplex algorithm
  tm <- try(
    nloptr(x0 = unname(x0_trans), 
           eval_f = ofun, 
           opts = list(algorithm = "NLOPT_LN_SBPLX", 
                       maxtime = 1 * 60, 
                       xtol_rel = 0.0,
                       print_level = 0))
  )
  
  # Display results
  print(tm$solution) # Parameters on the estimation scale
  print(-tm$objective) # Likelihood
}

#######################################################################################################
# End
#######################################################################################################