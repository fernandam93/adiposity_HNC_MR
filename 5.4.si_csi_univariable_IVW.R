#######################################################################
#               HEADSpAcE MENDELIAN RANDOMIZATION ANALYSES            #
#######################################################################
# R version 4.4.0 (2024-04-24)
# Last modified: 22 June 2025

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

#---------------------------------------------------------------------#
#                           Read exposure data                        #----
#---------------------------------------------------------------------#
#SI and CSI
si <- fread("results/clumped_exposure_files/clumped_si_incUKB.csv")
si <- si %>% filter(., pval.exposure <5e-8) #458 snps

csi <- fread("results/clumped_exposure_files/clumped_csi.csv")
csi <- csi %>% filter(., pval.exposure <5e-8) #458 snps

#---------------------------------------------------------------------#
#                   Prepare outcome data for analysis                #----
#---------------------------------------------------------------------#

#2. Extract SNP list for exposures from HEADSpAcE GWAS summary data 

source("scripts/alternative_functions.R")

#2.a.2 if outcome is in MR-BASE, this is the best alternative
#BMI (noUKB)
si_dat <- outcome_harmonisation_func_mrbase('ieu-b-5123',"HNC_noukb", si, "si_forFMB")
csi_dat <- outcome_harmonisation_func_mrbase('ieu-b-5123',"HNC_noukb", csi, "csi_forFMB")


#bind datasets
dat <- smartbind(si_dat, csi_dat) 


#---------------------------------------------------------------------#
#                   Mendelian randomization analysis                  #----
#---------------------------------------------------------------------#

mr_results <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", 
                                      "mr_weighted_median", "mr_weighted_mode"))
mr_results

mr_results[mr_results$exposure=="csi_forFMB",]$b <- mr_results[mr_results$exposure=="csi_forFMB",]$b*0.6940093
mr_results[mr_results$exposure=="csi_forFMB",]$se <- mr_results[mr_results$exposure=="csi_forFMB",]$se*0.6940093

#convert MR results to odds ratios
or_results <- generate_odds_ratios(mr_results)
or_results

results <- or_results %>%
  select(., exposure, outcome, method, nsnp, b, se, pval, or, or_lci95, or_uci95)

#save results
write.table(results, "results/MR_results/smoking_FMB_hnc_noukb_results.txt", row.names = F) 
