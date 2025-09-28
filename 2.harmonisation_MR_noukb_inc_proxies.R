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

#---------------------------------------------------------------------#
#                           Read exposure data                        #----
#---------------------------------------------------------------------#
#BMI Pulit
bmi_dat <- fread("results/clumped_exposure_files/clumped_BMI.csv")
bmi_dat2 <- bmi_dat %>% filter(., pval.exposure <5e-9) #458 snps

#WHR Pulit
whr_dat <- fread("results/clumped_exposure_files/clumped_WHR.csv")
whr_dat2 <- whr_dat %>% filter(., pval.exposure <5e-9) #283 snps

#WC Shungin (these were used originally, but analyses were updated to include a larger GWAS for waist circumference later)
#wc_dat <- fread("results/clumped_exposure_files/clumped_WC.csv") #41 snps

#WC UKB
wc_dat <- fread("results/clumped_exposure_files/clumped_WC2.csv") %>% mutate(, exposure = "WC") #375 snps

#---------------------------------------------------------------------#
#                   Prepare outcome data for analysis                #----
#---------------------------------------------------------------------#

#2. Extract SNP list for exposures from HEADSpAcE GWAS summary data 

source("scripts/alternative_functions.R")

#2.a.2 if outcome is in MR-BASE, this is the best alternative
#BMI (noUKB)
BMI_HNC_dat2 <- outcome_harmonisation_func_mrbase('ieu-b-5123',"HNC_noukb", bmi_dat2, "BMI_strict")
BMI_HPC_dat2 <- outcome_harmonisation_func_mrbase('ieu-b-5124',"HPC_noukb", bmi_dat2, "BMI_strict")
BMI_LA_dat2 <- outcome_harmonisation_func_mrbase('ieu-b-5125',"LA_noukb", bmi_dat2, "BMI_strict")
BMI_OC_dat2 <- outcome_harmonisation_func_mrbase('ieu-b-5126',"OC_noukb", bmi_dat2, "BMI_strict")
BMI_OPC_NEG_dat2 <- outcome_harmonisation_func_mrbase('ieu-b-5127',"OPC_NEG_noukb", bmi_dat2, "BMI_strict") 
BMI_OPC_POS_dat2 <- outcome_harmonisation_func_mrbase('ieu-b-5128',"OPC_POS_noukb", bmi_dat2, "BMI_strict")

#whr (noUKB)
WHR_HNC_dat2 <- outcome_harmonisation_func_mrbase('ieu-b-5123',"HNC_noukb", whr_dat2, "WHR_strict")
WHR_HPC_dat2 <- outcome_harmonisation_func_mrbase('ieu-b-5124',"HPC_noukb", whr_dat2, "WHR_strict")
WHR_LA_dat2 <- outcome_harmonisation_func_mrbase('ieu-b-5125',"LA_noukb", whr_dat2, "WHR_strict")
WHR_OC_dat2 <- outcome_harmonisation_func_mrbase('ieu-b-5126',"OC_noukb", whr_dat2, "WHR_strict")
WHR_OPC_NEG_dat2 <- outcome_harmonisation_func_mrbase('ieu-b-5127',"OPC_NEG_noukb", whr_dat2, "WHR_strict") 
WHR_OPC_POS_dat2 <- outcome_harmonisation_func_mrbase('ieu-b-5128',"OPC_POS_noukb", whr_dat2, "WHR_strict")

#wc (noUKB)
WC_HNC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5123',"HNC_noukb", wc_dat, "WC")
WC_HPC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5124',"HPC_noukb", wc_dat, "WC")
WC_LA_dat <- outcome_harmonisation_func_mrbase('ieu-b-5125',"LA_noukb", wc_dat, "WC")
WC_OC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5126',"OC_noukb", wc_dat, "WC")
WC_OPC_NEG_dat <- outcome_harmonisation_func_mrbase('ieu-b-5127',"OPC_NEG_noukb", wc_dat, "WC") 
WC_OPC_POS_dat <- outcome_harmonisation_func_mrbase('ieu-b-5128',"OPC_POS_noukb", wc_dat, "WC")


#bind datasets
dat <- smartbind(BMI_HNC_dat2, BMI_HPC_dat2, BMI_LA_dat2, BMI_OC_dat2, BMI_OPC_NEG_dat2, BMI_OPC_POS_dat2,
                 WHR_HNC_dat2, WHR_HPC_dat2, WHR_LA_dat2, WHR_OC_dat2, WHR_OPC_NEG_dat2, WHR_OPC_POS_dat2,
                 WC_HNC_dat, WC_HPC_dat, WC_LA_dat, WC_OC_dat, WC_OPC_NEG_dat, WC_OPC_POS_dat) 


#---------------------------------------------------------------------#
#                   Mendelian randomization analysis                  #----
#---------------------------------------------------------------------#

mr_results <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", 
                                      "mr_weighted_median", "mr_weighted_mode"))
mr_results

#convert MR results to odds ratios
or_results <- generate_odds_ratios(mr_results)
or_results

results <- or_results %>%
                 select(., exposure, outcome, method, nsnp, b, se, pval, or, or_lci95, or_uci95)

#save results
write.table(results[results$exposure=="BMI_strict",], "results/MR_results/BMI_strict_hnc_noukb_results.txt", row.names = F) 
write.table(results[results$exposure=="WHR_strict",], "results/MR_results/WHR_strict_hnc_noukb_results.txt", row.names = F) 
write.table(results[results$exposure=="WC",], "results/MR_results/WC2_hnc_noukb_results.txt", row.names = F) 


