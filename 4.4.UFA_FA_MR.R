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
#Favourable and unfavourable adiposity
fa_dat <- fread("results/clumped_exposure_files/clumped_FA.csv")
ufa_dat <- fread("results/clumped_exposure_files/clumped_UFA.csv")

#body fat percentage in UKB
exp_dat <- extract_instruments(
  c("ukb-b-8909"),
  p1 = 5e-08,
  clump = TRUE,
  p2 = 5e-08,
  r2 = 0.001,
  kb = 10000,
  opengwas_jwt = ieugwasr::get_opengwas_jwt(),
  force_server = FALSE
)

bodyfat_dat <- clump_data(exp_dat, clump_p1 = 5e-8, clump_p2 = 5e-8)


#---------------------------------------------------------------------#
#                     Prepare outcomes for analysis                   #----
#---------------------------------------------------------------------#

source("scripts/alternative_functions.R")

#2.a.2 if outcome is in MR-BASE, this is the best alternative
#FA
FA_HNC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5123',"HNC_noukb", fa_dat, "FA")
FA_HPC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5124',"HPC_noukb", fa_dat, "FA")
FA_LA_dat <- outcome_harmonisation_func_mrbase('ieu-b-5125',"LA_noukb", fa_dat, "FA")
FA_OC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5126',"OC_noukb", fa_dat, "FA")
FA_OPC_NEG_dat <- outcome_harmonisation_func_mrbase('ieu-b-5127',"OPC_NEG_noukb", fa_dat, "FA") 
FA_OPC_POS_dat <- outcome_harmonisation_func_mrbase('ieu-b-5128',"OPC_POS_noukb", fa_dat, "FA")


#UFA
UFA_HNC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5123',"HNC_noukb", ufa_dat, "UFA")
UFA_HPC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5124',"HPC_noukb", ufa_dat, "UFA")
UFA_LA_dat <- outcome_harmonisation_func_mrbase('ieu-b-5125',"LA_noukb", ufa_dat, "UFA")
UFA_OC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5126',"OC_noukb", ufa_dat, "UFA")
UFA_OPC_NEG_dat <- outcome_harmonisation_func_mrbase('ieu-b-5127',"OPC_NEG_noukb", ufa_dat, "UFA") 
UFA_OPC_POS_dat <- outcome_harmonisation_func_mrbase('ieu-b-5128',"OPC_POS_noukb", ufa_dat, "UFA")


#body fat
bodyfat_HNC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5123',"HNC_noukb", bodyfat_dat, "Body fat percentage")
bodyfat_HPC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5124',"HPC_noukb", bodyfat_dat, "Body fat percentage")
bodyfat_LA_dat <- outcome_harmonisation_func_mrbase('ieu-b-5125',"LA_noukb", bodyfat_dat, "Body fat percentage")
bodyfat_OC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5126',"OC_noukb", bodyfat_dat, "Body fat percentage")
bodyfat_OPC_NEG_dat <- outcome_harmonisation_func_mrbase('ieu-b-5127',"OPC_NEG_noukb", bodyfat_dat, "Body fat percentage") 
bodyfat_OPC_POS_dat <- outcome_harmonisation_func_mrbase('ieu-b-5128',"OPC_POS_noukb", bodyfat_dat, "Body fat percentage")


dat <- smartbind(FA_HNC_dat, FA_HPC_dat, FA_LA_dat, FA_OC_dat, FA_OPC_NEG_dat, FA_OPC_POS_dat,
                 UFA_HNC_dat, UFA_HPC_dat, UFA_LA_dat, UFA_OC_dat, UFA_OPC_NEG_dat, UFA_OPC_POS_dat,
                 bodyfat_HNC_dat, bodyfat_HPC_dat, bodyfat_LA_dat, bodyfat_OC_dat, bodyfat_OPC_NEG_dat, bodyfat_OPC_POS_dat)


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

write.table(results, "results/MR_results/secondary_adiposity/bodyfat_results.txt", row.names = F) 
#mr_results <- read.table("results/MR_results/secondary_adiposity/bodyfat_results.txt", header = T)
#---------------------------------------------------------------------#
#                         Prepare results for plots                   #----
#---------------------------------------------------------------------#
mr_results$outcome[mr_results$outcome=="OC_noukb"] <- "Oral cavity cancer"
mr_results$outcome[mr_results$outcome=="HPC_noukb"] <- "Hypopharyngeal cancer"
mr_results$outcome[mr_results$outcome=="OPC_POS_noukb"] <- "HPV positive oropharyngeal cancer"
mr_results$outcome[mr_results$outcome=="OPC_NEG_noukb"] <- "HPV negative oropharyngeal cancer"
mr_results$outcome[mr_results$outcome=="LA_noukb"] <- "Laryngeal cancer"
mr_results$outcome[mr_results$outcome=="HNC_noukb"] <- "Head and neck cancer"

mr_results$exposure[mr_results$exposure=="FA"] <- "Favourable adiposity"
mr_results$exposure[mr_results$exposure=="UFA"] <- "Unfavourable adiposity"

cancers <- c('Head and neck cancer', 'Oral cavity cancer', 'Laryngeal cancer', 
             'Hypopharyngeal cancer', 'HPV positive oropharyngeal cancer', 'HPV negative oropharyngeal cancer')

mr_results <- split_exposure(mr_results)

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
    xlim= c(0.01,30)
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

combine2_func <- function(A, B){
  plot.new()
  par(mar=c(1,1,1,1), mgp=c(3,1,0))
  x <- ggarrange(A, B, labels=c('A', 'B'),
                 ncol = 2, nrow = 1, common.legend = T, hjust = -3, legend = "bottom")
  print(x)
}

#---------------------------------------------------------------------#
#                            Forest plots                             #----
#---------------------------------------------------------------------#
p_fa <- myforestplot(mr_results[mr_results$exposure=="Favourable adiposity",], "", "Odds ratio (95% CI) per 1-SD higher favourable adiposity")
p_ufa <- myforestplot(mr_results[mr_results$exposure=="Unfavourable adiposity",], "", "Odds ratio (95% CI) per 1-SD higher unfavourable adiposity")
p_bodyfat <- myforestplot(mr_results[mr_results$exposure=="Body fat percentage",], "", "Odds ratio (95% CI) per 1-SD higher body fat percentage")


new <- combine2_func(p_fa, p_ufa)

#---------------------------------------------------------------------#
#                  Save results and plots                             #----
#---------------------------------------------------------------------#

#save plot function
save_forest_func <- function(file_name, plot_name)
{
  png(file_name, res=330, height=2000, width=3000)
  print(plot_name)
  dev.off()
}
#save
save_forest_func('results/plots/secondary_adiposity/bodyfat_plot.png', p_bodyfat)


#save plot function2
save_func <- function(file_name, plot_name)
{
  png(file_name, res=330, height=2500, width=5000)
  print(plot_name)
  dev.off()
}
#save
save_func('results/plots/secondary_adiposity/FA_UFA_plot.png', new)

