#######################################################################################################
# This script installs the packages needed to run the other scripts in the folder
# It uses the 'checkpoint' package to set the package versions to those that were used for analysis
# The checkpoint date is set to 1 October 2015, which corresponds to pomp version 1.2.1.1
#######################################################################################################

rm(list = ls())
ipkg <- rownames(installed.packages()) # List installed packages

# Check if the 'checkpoint' package is already installed and install it if necessary
if(!("checkpoint" %in% ipkg)) install.packages("checkpoint")

library(checkpoint)
checkpoint(snapshotDate = "2015-10-01", scanForPackages = T)