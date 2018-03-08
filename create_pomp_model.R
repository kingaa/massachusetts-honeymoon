#######################################################################################################
# This function creates the pomp object for the age-structured VSEIR2 model
#######################################################################################################
create_pomp_model <- function(data, covars, nages_mod, nbasis, time_start_data, time_start_sim, dt) {
  # Args:
  # data: data frame of observations, with time and no of cases in each age group (7 groups)
  # covars: data frame of covariates: time, yearly no of births, and age-specific migration rate
  # nages_mod: no of age groups in the model 
  # nseas: no of periodic spline bases to model seasonality in children and teenagers
  # time_start_data: time at which the data start
  # time_start_sim: time at which the simulations start; should be before the start time for data,
  # and after the start date of the table of covariates
  # dt: time step (in years) for stochastic simulations
  # Returns: pomp object
  
  # Create the shared library from the compiled code ------------------------
  model <- "model_equations"
  modelfile <- paste0(model, ".c")
  solib <- paste0(model, .Platform$dynlib.ext)
  
  ## sets PKG_CFLAGS to include path to pomp header file
  ## [This won't work on windoze.]
  cflags <- paste0("PKG_CFLAGS=\"",
                   Sys.getenv("PKG_CFLAGS"),
                   " -I", system.file("include", package = "pomp"),"\"")
  
  ## uses 'system2' instead of 'system'
  rv <- system2(
    command = R.home("bin/R"),
    args = c("CMD","SHLIB","-o", solib, modelfile),
    env = cflags)
  if(rv != 0) stop("cannot compile shared-object library ", solib)
  dyn.load(solib)
  
  # Create spline bases
  # Each basis had period 1y, and mean 1 / nseas over one year
  seas <- periodic.bspline.basis(x = covars$time, nbasis = nbasis, degree = 3, period = 1, names = "seas%d")
  covars <- data.frame(covars, seas)
  
  # Rename data
  colnames(data)[-1] <- paste0("reports-", 1:(ncol(data) - 1))
  
  # Shift time for data and covariates
  data$time <- data$time - time_start_data
  covars$time <- covars$time - time_start_data
  
  po <- pomp(
    data = data,
    times = "time",
    t0 = time_start_sim,
    obsnames = colnames(data)[-1],
    covar = covars, 
    tcovar = "time",
    covarnames = c("births", "mu1", "seas1"),
    statenames = paste0(c("V", "S1", "E1", "I1", "S2", "E2", "I2", "R", "C1", "C2"), "-1"),
    zeronames = c(paste0("C1-", 1:nages_mod), paste0("C2-", 1:nages_mod)),
    paramnames = c("iota", "q1-1", "q2-1", "omegaC1", "omegaT1",
                   "CR1", "theta", "beta_sd", "sigma", "gamma", 
                   "epsilonV", "epsilonI", "alphaV", "alphaI", "rho1-1", 
                   "rho2", "tau", "epsA", "v1", "v2",
                   "t0", "t1", "t2", "delta1", "N1",
                   "model_rho", "model_q", "model_vac",
                   as.character(sapply(paste0(c("V", "S1", "E1", "I1", "S2", "E2", "I2", "R"), "-"), 
                                       paste0, 1:nages_mod, "-0"))),
    PACKAGE = "model_equations",
    dmeasure = "dObs",
    rmeasure = "rObs",
    rprocess = euler.sim(step.fun = "rSim", delta.t = dt),
    skeleton = "skel_continuous",
    skeleton.type = "vectorfield",
    fromEstimationScale = "par_trans",
    toEstimationScale = "par_untrans",
    nseas = nbasis,
    initializer = function(params, t0, ...) {
      all.states <- paste0(c("V", "S1", "E1", "I1", "S2", "E2", "I2", "R", "C1", "C2"), "-") # All state variables
      comp.states <- paste0(c("V", "S1", "E1", "I1", "S2", "E2", "I2", "R"), "-") # All compartments
      all.state.names <- sapply(all.states, paste0, 1:nages_mod) %>% as.character(.)
      comp.names <- sapply(comp.states, paste0, 1:nages_mod) %>% as.character(.)
      comp.ic.names <- sapply(comp.states, paste0, 1:nages_mod, "-0") %>% as.character(.)
      x0 <- setNames(numeric(length(all.state.names)), all.state.names) 
      frac <- params[comp.ic.names] # Initial fractions in each compartment
      pop <- rep(params[paste0("N", 1:nages_mod)], length(comp.states)) # Population sizes in each age group, repeated nages times
      x0[comp.names] <- round(frac * pop)
      x0
    }
  )
  return(po)
}