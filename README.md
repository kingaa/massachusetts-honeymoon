# massachusetts-honeymoon
Supporting material for the 2018 *Science Translational Medicine* paper, "The impact of past vaccination coverage and immunity on pertussis resurgence", by M. Domenech de Cell&egrave;s, F. M. G. Magpantay, A. A. King, and P. Rohani.


This repository contains saved the **R** objects (`.rds`) and **R** scripts (`.R`) needed to run the estimations for the deterministic and stochastic models.
[Download](https://github.com/kingaa/massachusetts-honeymoon/archive/master.zip) and unpack all the files first.
Next, run the script `install_packages.R`;
it installs the "checkpoint" package (if needed) and the other required packages at their version on 1 October 2015.
(In particular, this ensures that the **pomp** package version is 1.2.1.1)

The main scripts are then:

1. `run_traj_match.R`: run trajectory matching for the deterministic models.
   This script example runs a single estimation from one parameter set starting value, with maximum execution time of 1 min.
2. `run_mifs.R`: run the maximum iterated filtering algorithm to estimate the parameters of the stochastic models.
   This script example runs a single estimation (with 2 MIF iterations and 2e2 particles), followed by 2 evaluations of the likelihood using particle filtering (with 2e2 particles).

Other items:

1. **R** objects:
     - `covars.rds`: data frame containing the covariates used in the model (birth and age-specific migration rates)
     - `data.rds`: data frame that contains the data, i.e., monthly age-specific case reports in Massachusetts during 1990–2005
     - `mle_waning_model_deterministic.rds`: named vector containing the model parameters (fixed+estimated MLE) of the deterministic model
     - `mle_waning_model_stochastic.rds`: named vector containing the model parameters (fixed+estimated MLE) of the stochastic model
2. Other files:
     - `model_equations.c`: C code implementing the observation and the process (deterministic and stochastic variants) models
     - `create_pomp_model.R`: **R** script implementing a function that creates the pomp object. 
