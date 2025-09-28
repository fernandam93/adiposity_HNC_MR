#######################################################################
#               HEADSpAcE MENDELIAN RANDOMIZATION ANALYSES            #
#######################################################################
# R version 4.4.0 (2024-04-24)
# Last modified: 17 Jan 2025

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("YOUR_WD") 

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot","gtools", "scales", "meta", "metafor")

# or
# install.packages("devtools")
# install.packages("tidyverse")
# devtools::install_github("NightingaleHealth/ggforestplot")
# devtools::install_github("r-gregmisc/gtools")
#---------------------------------------------------------------------#
#                             Read MR results                         #----
#---------------------------------------------------------------------#

results1_noukb <- read.table("results/MR_results/WHR_strict_hnc_noukb_results.txt", header = T) 
results2_noukb <- read.table("results/MR_results/BMI_strict_hnc_noukb_results.txt", header = T) 
results3_noukb <- read.table("results/MR_results/WC2_hnc_noukb_results.txt", header = T) 


results <- rbind(results1_noukb, results2_noukb, results3_noukb)


# Function to run a fixed effect meta-analysis 
meta_func <- function(method_varname, exp_varname, out1, out2="", out3="", out4="", 
                      out5="", out6="", out7="", out8="", out9="", out10="", 
                      out11="", out12="", out13="", out14="", out15="", 
                      out16="", out17="", out18="", out19="", out20="", out21="")
{
  input <- results[which(results$method==method_varname),]
  input <- input[which(input$exposure==exp_varname),]
  input <- input[which(input$outcome==out1 | input$outcome==out2 | input$outcome==out3 | 
                         input$outcome==out4 | input$outcome==out5 | input$outcome==out6 | 
                         input$outcome==out7 | input$outcome==out8 | input$outcome==out9 | 
                         input$outcome==out10 | input$outcome==out11 | input$outcome==out12 |
                         input$outcome==out13 | input$outcome==out14 | input$outcome==out15 | 
                         input$outcome==out16 | input$outcome==out17 | input$outcome==out18 | 
                         input$outcome==out19 | input$outcome==out20 | input$outcome==out21),]
  #meta-analysis
  a <- metagen(TE = b, seTE = se, data = input, 
               studlab = paste(outcome), sm = "OR",
               hakn = FALSE, byvar = c(exposure),
               method.tau="DL", comb.fixed = T, comb.random = F) 
  print(a)
  
  #extract values from meta output
  TE.tibble <- as_tibble(a$TE.fixed.w)
  se.tibble <- as_tibble(a$seTE.fixed.w)
  p.tibble <- as_tibble(a$pval.fixed.w)
  bylevs.tibble <- as_tibble(a$bylevs)
  #combine tibbles and change column names
  tibble <- cbind(TE.tibble, se.tibble, p.tibble, bylevs.tibble)
  colnames(tibble) <- c("b", "se", "pval", "exposure")
  #add columns for exposure and method
  tibble$outcomes <- paste(out1, out2, out3, out4, out5, sep = ", ")
  tibble$het_test <- a$pval.Q
  tibble$method <- method_varname
  tibble
  
}

#IVW
IVW_bmi_noukb <- meta_func("Inverse variance weighted", "BMI_strict", "OC_noukb", "LA_noukb", "HPC_noukb", 
                         "OPC_NEG_noukb", "OPC_POS_noukb")
IVW_whr_noukb <- meta_func("Inverse variance weighted", "WHR_strict", "OC_noukb", "LA_noukb", "HPC_noukb", 
                     "OPC_NEG_noukb", "OPC_POS_noukb")
IVW_wc_noukb <- meta_func("Inverse variance weighted", "WC", "OC_noukb", "LA_noukb", "HPC_noukb", 
                     "OPC_NEG_noukb", "OPC_POS_noukb")

