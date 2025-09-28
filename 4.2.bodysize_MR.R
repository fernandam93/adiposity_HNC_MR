#######################################################################
#               HEADSpAcE MENDELIAN RANDOMIZATION ANALYSES            #
#######################################################################
# R version 4.4.0 (2024-04-24)
# Last modified: 5 August 2024

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
#Richardson body size
adult_bmi <- extract_instruments(
            c("ieu-b-5118"),
            p1 = 5e-08,
            clump = TRUE,
            p2 = 5e-08,
            r2 = 0.001,
            kb = 10000,
            opengwas_jwt = ieugwasr::get_opengwas_jwt(),
            force_server = FALSE
          )

adult_bmi <- clump_data(adult_bmi, clump_p1 = 5e-8, clump_p2 = 5e-8)

child_bmi <- extract_instruments(
  c("ieu-b-5107"),
  p1 = 5e-08,
  clump = TRUE,
  p2 = 5e-08,
  r2 = 0.001,
  kb = 10000,
  opengwas_jwt = ieugwasr::get_opengwas_jwt(),
  force_server = FALSE
)

child_bmi <- clump_data(child_bmi, clump_p1 = 5e-8, clump_p2 = 5e-8)

#---------------------------------------------------------------------#
#                     Prepare outcomes for analysis                   #----
#---------------------------------------------------------------------#

source("scripts/alternative_functions.R")

#2.a.2 if outcome is in MR-BASE, this is the best alternative
#body_size
adult_body_size_HNC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5123',"HNC_noukb", adult_bmi, "adult BMI")
adult_body_size_HPC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5124',"HPC_noukb", adult_bmi, "adult BMI")
adult_body_size_LA_dat <- outcome_harmonisation_func_mrbase('ieu-b-5125',"LA_noukb", adult_bmi, "adult BMI")
adult_body_size_OC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5126',"OC_noukb", adult_bmi, "adult BMI")
adult_body_size_OPC_NEG_dat <- outcome_harmonisation_func_mrbase('ieu-b-5127',"OPC_NEG_noukb", adult_bmi, "adult BMI") 
adult_body_size_OPC_POS_dat <- outcome_harmonisation_func_mrbase('ieu-b-5128',"OPC_POS_noukb", adult_bmi, "adult BMI")

child_body_size_HNC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5123',"HNC_noukb", child_bmi, "child BMI")
child_body_size_HPC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5124',"HPC_noukb", child_bmi, "child BMI")
child_body_size_LA_dat <- outcome_harmonisation_func_mrbase('ieu-b-5125',"LA_noukb", child_bmi, "child BMI")
child_body_size_OC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5126',"OC_noukb", child_bmi, "child BMI")
child_body_size_OPC_NEG_dat <- outcome_harmonisation_func_mrbase('ieu-b-5127',"OPC_NEG_noukb", child_bmi, "child BMI") 
child_body_size_OPC_POS_dat <- outcome_harmonisation_func_mrbase('ieu-b-5128',"OPC_POS_noukb", child_bmi, "child BMI")

dat <- smartbind(adult_body_size_HNC_dat, adult_body_size_HPC_dat, adult_body_size_LA_dat, 
                 adult_body_size_OC_dat, adult_body_size_OPC_NEG_dat, adult_body_size_OPC_POS_dat,
                 child_body_size_HNC_dat, child_body_size_HPC_dat, child_body_size_LA_dat, 
                 child_body_size_OC_dat, child_body_size_OPC_NEG_dat, child_body_size_OPC_POS_dat)


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

write.table(results, "results/MR_results/secondary_adiposity/bodysize_univar_results.txt", row.names = F) 
#mr_results <- read.table("results/MR_results/secondary_adiposity/bodysize_univar_results.txt", header = T)

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
    xlim= c(0,5)
  ) +
    ggforce::facet_col(
      facets = ~group,
      scales = "free_y",
      space = "free"
    )
  colours_BP <- c("#999999", "#E69F00", "#56B4E9", "#009E73")
  x <- x + scale_color_manual(values=colours_BP)
  x <- x + scale_x_continuous(breaks = c(0.3, 1, 3, 5))
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
p_bodysize <- myforestplot(mr_results[mr_results$exposure=="adult BMI",], "", "Odds ratio (95% CI) per change in adult body size category")
p_bodysize_child <- myforestplot(mr_results[mr_results$exposure=="child BMI",], "", "Odds ratio (95% CI) per change in child body size category")

new <- combine2_func(p_bodysize_child, p_bodysize)

#---------------------------------------------------------------------#
#                  Save results and plots                             #----
#---------------------------------------------------------------------#

#save plot
save_func <- function(file_name, plot_name)
{
  png(file_name, res=330, height=2500, width=5000)
  print(plot_name)
  dev.off()
}

#save

save_func('results/plots/secondary_adiposity/bodysize_plot.png', new)


