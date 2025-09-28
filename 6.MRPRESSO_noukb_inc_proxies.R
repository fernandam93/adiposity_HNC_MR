#######################################################################
#               HEADSpAcE MENDELIAN RANDOMIZATION ANALYSES            #
#######################################################################
# R version 4.4.0 (2024-04-24)
# Last modified: 17 January 2025

#script to harmonise exposure and outcome datasets in MR-Base (excluding UKB from outcome)
#requires scripts/alternative_functions.R

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment 
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("YOUR_WD") 

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot","gtools", "LDlinkR", "simex", "ieugwasr")
#pacman::p_load_gh("MRCIEU/MRInstruments", "MRCIEU/TwoSampleMR") #use this if you haven't installed these packages before

#alternative way of installing and loading MRCIEU packages from GitHub
# install.packages("devtools")
# devtools::install_github("MRCIEU/TwoSampleMR")
# devtools::install_github("MRCIEU/MRInstruments", force = TRUE)
# library(TwoSampleMR)
# library(MRInstruments)
#devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)

#---------------------------------------------------------------------#
#                           Read harmonised data                        #----
#---------------------------------------------------------------------#

WHR_HNCnoukb <- read.csv("results/harmonised_data_files/harmonised_WHR_strict_HNC_noukb.csv", header = T) 
BMI_HNCnoukb <- read.csv("results/harmonised_data_files/harmonised_BMI_strict_HNC_noukb.csv", header = T) 
WC_HNCnoukb <- read.csv("results/harmonised_data_files/harmonised_WC_HNC_noukb.csv", header = T)

#---------------------------------------------------------------------#
#                             MR-PRESSO analysis                       #----
#---------------------------------------------------------------------#

mrpresso_results_WHR_HNC <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = WHR_HNCnoukb, NbDistribution = 6000,  SignifThreshold = 0.05, seed = 123123)
mrpresso_results_BMI_HNC <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = BMI_HNCnoukb, NbDistribution = 10000,  SignifThreshold = 0.05, seed = 123123)
mrpresso_results_WC_HNC <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = WC_HNCnoukb, NbDistribution = 8000,  SignifThreshold = 0.05, seed = 123123)


capture.output(mrpresso_results_WHR_HNC, file = "results/MR_results/mr-presso/mrpresso_WHR_HNC.txt")
capture.output(mrpresso_results_BMI_HNC, file = "results/MR_results/mr-presso/mrpresso_BMI_HNC.txt")
capture.output(mrpresso_results_WC_HNC, file = "results/MR_results/mr-presso/mrpresso_WC2_HNC.txt")







