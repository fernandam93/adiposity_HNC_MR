#######################################################################
#               HEADSpAcE MENDELIAN RANDOMIZATION ANALYSES            #
#######################################################################
# R version 4.4.0 (2024-04-24)
# Last modified: 16 June 2024

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
#GIANT body shape PCs
pc1_dat <- fread("results/clumped_exposure_files/clumped_PC1.csv")
pc2_dat <- fread("results/clumped_exposure_files/clumped_PC2.csv")
pc3_dat <- fread("results/clumped_exposure_files/clumped_PC3.csv")
pc4_dat <- fread("results/clumped_exposure_files/clumped_PC4.csv")


#---------------------------------------------------------------------#
#                     Prepare outcomes for analysis                   #----
#---------------------------------------------------------------------#

source("scripts/alternative_functions.R")

#2.a.2 if outcome is in MR-BASE, this is the best alternative
#body_shape PC1
PC1_HNC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5123',"HNC_noukb", pc1_dat, "PC1")
PC1_HPC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5124',"HPC_noukb", pc1_dat, "PC1")
PC1_LA_dat <- outcome_harmonisation_func_mrbase('ieu-b-5125',"LA_noukb", pc1_dat, "PC1")
PC1_OC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5126',"OC_noukb", pc1_dat, "PC1")
PC1_OPC_NEG_dat <- outcome_harmonisation_func_mrbase('ieu-b-5127',"OPC_NEG_noukb", pc1_dat, "PC1") 
PC1_OPC_POS_dat <- outcome_harmonisation_func_mrbase('ieu-b-5128',"OPC_POS_noukb", pc1_dat, "PC1")


#body_shape PC2
PC2_HNC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5123',"HNC_noukb", pc2_dat, "PC2")
PC2_HPC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5124',"HPC_noukb", pc2_dat, "PC2")
PC2_LA_dat <- outcome_harmonisation_func_mrbase('ieu-b-5125',"LA_noukb", pc2_dat, "PC2")
PC2_OC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5126',"OC_noukb", pc2_dat, "PC2")
PC2_OPC_NEG_dat <- outcome_harmonisation_func_mrbase('ieu-b-5127',"OPC_NEG_noukb", pc2_dat, "PC2") 
PC2_OPC_POS_dat <- outcome_harmonisation_func_mrbase('ieu-b-5128',"OPC_POS_noukb", pc2_dat, "PC2")


#body_shape PC3
PC3_HNC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5123',"HNC_noukb", pc3_dat, "PC3")
PC3_HPC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5124',"HPC_noukb", pc3_dat, "PC3")
PC3_LA_dat <- outcome_harmonisation_func_mrbase('ieu-b-5125',"LA_noukb", pc3_dat, "PC3")
PC3_OC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5126',"OC_noukb", pc3_dat, "PC3")
PC3_OPC_NEG_dat <- outcome_harmonisation_func_mrbase('ieu-b-5127',"OPC_NEG_noukb", pc3_dat, "PC3") 
PC3_OPC_POS_dat <- outcome_harmonisation_func_mrbase('ieu-b-5128',"OPC_POS_noukb", pc3_dat, "PC3")


#body_shape PC4
PC4_HNC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5123',"HNC_noukb", pc4_dat, "PC4")
PC4_HPC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5124',"HPC_noukb", pc4_dat, "PC4")
PC4_LA_dat <- outcome_harmonisation_func_mrbase('ieu-b-5125',"LA_noukb", pc4_dat, "PC4")
PC4_OC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5126',"OC_noukb", pc4_dat, "PC4")
PC4_OPC_NEG_dat <- outcome_harmonisation_func_mrbase('ieu-b-5127',"OPC_NEG_noukb", pc4_dat, "PC4") 
PC4_OPC_POS_dat <- outcome_harmonisation_func_mrbase('ieu-b-5128',"OPC_POS_noukb", pc4_dat, "PC4")



dat <- smartbind(PC1_HNC_dat, PC1_HPC_dat, PC1_LA_dat, PC1_OC_dat, PC1_OPC_NEG_dat, PC1_OPC_POS_dat,
                 PC2_HNC_dat, PC2_HPC_dat, PC2_LA_dat, PC2_OC_dat, PC2_OPC_NEG_dat, PC2_OPC_POS_dat,
                 PC3_HNC_dat, PC3_HPC_dat, PC3_LA_dat, PC3_OC_dat, PC3_OPC_NEG_dat, PC3_OPC_POS_dat,
                 PC4_HNC_dat, PC4_HPC_dat, PC4_LA_dat, PC4_OC_dat, PC4_OPC_NEG_dat, PC4_OPC_POS_dat)


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

write.table(results, "results/MR_results/secondary_adiposity/bodyshape_results.txt", row.names = F) 

#---------------------------------------------------------------------#
#                         Prepare results for plots                   #----
#---------------------------------------------------------------------#
mr_results$outcome[mr_results$outcome=="OC_noukb"] <- "Oral cavity cancer"
mr_results$outcome[mr_results$outcome=="HPC_noukb"] <- "Hypopharyngeal cancer"
mr_results$outcome[mr_results$outcome=="OPC_POS_noukb"] <- "HPV positive oropharyngeal cancer"
mr_results$outcome[mr_results$outcome=="OPC_NEG_noukb"] <- "HPV negative oropharyngeal cancer"
mr_results$outcome[mr_results$outcome=="LA_noukb"] <- "Laryngeal cancer"
mr_results$outcome[mr_results$outcome=="HNC_noukb"] <- "Head and neck cancer"


cancers <- c('Head and neck cancer', 'Oral cavity cancer', 'Laryngeal cancer', 
             'Hypopharyngeal cancer', 'HPV positive oropharyngeal cancer', 'HPV negative oropharyngeal cancer')

as.factor(mr_results$outcome)
mr_results$outcome <- factor(mr_results$outcome, levels = cancers)
mr_results <- arrange(mr_results, outcome)


as.factor(mr_results$method)
mr_results$method <- factor(mr_results$method, levels = c("Weighted mode", "Weighted median",
                                                                    "MR Egger", "Inverse variance weighted"))

mr_results$group[mr_results$outcome!="Head and neck cancer"] <- "Subsites"
mr_results$group[mr_results$outcome=="Head and neck cancer"] <- "Overall"


#---------------------------------------------------------------------#
#                            Forest function                          #----
#---------------------------------------------------------------------#
#plots
myforestplot <- function(df, exp_dataset, xlab)
{
  options(scipen = 999)
  x <- forestplot(
    df = df,
    estimate = b,
    se = se,
    pvalue = pval,
    name = outcome,
    logodds = T,
    colour = method,
    title = exp_dataset,
    xlab = xlab,
    #xlim= c(0.1,5)
  ) +
    ggforce::facet_col(
      facets = ~group,
      scales = "free_y",
      space = "free"
    )
  colours_BP <- c("#999999", "#E69F00", "#56B4E9", "#009E73")
  x <- x + scale_color_manual(values=colours_BP)
  #x <- x + scale_x_continuous(breaks = c(0.3, 1, 3, 5))
  print(x)
}

combine2_func <- function(A, B, C, D){
  plot.new()
  par(mar=c(1,1,1,1), mgp=c(3,1,0))
  x <- ggarrange(A, B, C, D, labels=c('A', 'B', 'C', 'D'),
                 ncol = 4, nrow = 1, common.legend = T, hjust = -3, legend = "bottom")
  print(x)
}

#---------------------------------------------------------------------#
#                            Forest plots                             #----
#---------------------------------------------------------------------#
p_pc1 <- myforestplot(mr_results[mr_results$exposure=="PC1",], "", "Odds ratio (95% CI) per 1-SD higher PC1")
p_pc2 <- myforestplot(mr_results[mr_results$exposure=="PC2",], "", "Odds ratio (95% CI) per 1-SD higher PC2")
p_pc3 <- myforestplot(mr_results[mr_results$exposure=="PC3",], "", "Odds ratio (95% CI) per 1-SD higher PC3")
p_pc4 <- myforestplot(mr_results[mr_results$exposure=="PC4",], "", "Odds ratio (95% CI) per 1-SD higher PC4")


new <- combine2_func(p_pc1, p_pc2, p_pc3, p_pc4)

#---------------------------------------------------------------------#
#                  Save results and plots                             #----
#---------------------------------------------------------------------#

#save plot
save_func <- function(file_name, plot_name)
{
  png(file_name, res=330, height=2500, width=8000)
  print(plot_name)
  dev.off()
}

#save

save_func('results/plots/secondary_adiposity/bodyshape_plot.png', new)


