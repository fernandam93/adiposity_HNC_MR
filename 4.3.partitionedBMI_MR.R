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
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot","gtools", "LDlinkR", "simex", "ieugwasr", "openxlsx")
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
#GIANT partitioned BMI
brain_dat <- read.xlsx("data/GIANT_Pulit2018/tissue_partitioned_BMI.xlsx") %>% filter(., Tissue=="brain" | Tissue=="both") %>% format_data(type = "exposure") 
adipose_dat <- read.xlsx("data/GIANT_Pulit2018/tissue_partitioned_BMI.xlsx") %>% filter(., Tissue=="adipose" | Tissue=="both") %>% format_data(type = "exposure") 

#---------------------------------------------------------------------#
#                     Prepare outcomes for analysis                   #----
#---------------------------------------------------------------------#

#2. Extract SNP list for exposures from HEADSpAcE GWAS summary data 

source("scripts/alternative_functions.R")

#2.a.2 if outcome is in MR-BASE, this is the best alternative
#adipose BMI
adi_HNC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5123',"HNC_noukb", adipose_dat, "adipose-BMI")
adi_HPC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5124',"HPC_noukb", adipose_dat, "adipose-BMI")
adi_LA_dat <- outcome_harmonisation_func_mrbase('ieu-b-5125',"LA_noukb", adipose_dat, "adipose-BMI")
adi_OC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5126',"OC_noukb", adipose_dat, "adipose-BMI")
adi_OPC_NEG_dat <- outcome_harmonisation_func_mrbase('ieu-b-5127',"OPC_NEG_noukb", adipose_dat, "adipose-BMI") 
adi_OPC_POS_dat <- outcome_harmonisation_func_mrbase('ieu-b-5128',"OPC_POS_noukb", adipose_dat, "adipose-BMI")


#brain BMI 
brain_HNC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5123',"HNC_noukb", brain_dat, "brain-BMI")
brain_HPC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5124',"HPC_noukb", brain_dat, "brain-BMI")
brain_LA_dat <- outcome_harmonisation_func_mrbase('ieu-b-5125',"LA_noukb", brain_dat, "brain-BMI")
brain_OC_dat <- outcome_harmonisation_func_mrbase('ieu-b-5126',"OC_noukb", brain_dat, "brain-BMI")
brain_OPC_NEG_dat <- outcome_harmonisation_func_mrbase('ieu-b-5127',"OPC_NEG_noukb", brain_dat, "brain-BMI") 
brain_OPC_POS_dat <- outcome_harmonisation_func_mrbase('ieu-b-5128',"OPC_POS_noukb", brain_dat, "brain-BMI")



dat <- smartbind(adi_HNC_dat, adi_HPC_dat, adi_LA_dat, adi_OC_dat, adi_OPC_NEG_dat, adi_OPC_POS_dat,
                 brain_HNC_dat, brain_HPC_dat, brain_LA_dat, brain_OC_dat, brain_OPC_NEG_dat, brain_OPC_POS_dat)


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

write.table(results, "results/MR_results/secondary_adiposity/tissuepartitioned_results.txt", row.names = F) 

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
    xlim= c(0.01,60)
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
p_adi <- myforestplot(mr_results[mr_results$exposure=="adipose-BMI",], "", "Odds ratio (95% CI) per 1-SD higher adipose BMI")
p_brain <- myforestplot(mr_results[mr_results$exposure=="brain-BMI",], "", "Odds ratio (95% CI) per 1-SD higher brain BMI")



new <- combine2_func(p_adi, p_brain)

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

save_func('results/plots/secondary_adiposity/tissuepartitioned_plot.png', new)


