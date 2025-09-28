#######################################################################
#               HEADSpAcE MENDELIAN RANDOMIZATION ANALYSES            #
#######################################################################
# R version 4.4.0 (2024-04-24)
# Last modified: 17 January 2025

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("YOUR_WD") 

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot","gtools", "scales")
 # or
# install.packages("devtools")
# install.packages("tidyverse")
# devtools::install_github("NightingaleHealth/ggforestplot")
# devtools::install_github("r-gregmisc/gtools")

#---------------------------------------------------------------------#
#                             Read MR results                         #----
#---------------------------------------------------------------------#

#read results
results_WHR_HNCnoukb <- read.table("results/MR_results/WHR_strict_hnc_noukb_results.txt", header = T) %>% mutate(., overlap = "non-overlapping")
results_BMI_HNCnoukb <- read.table("results/MR_results/BMI_strict_hnc_noukb_results.txt", header = T) %>% mutate(., overlap = "non-overlapping")
results_WC_HNCnoukb <- read.table("results/MR_results/WC2_hnc_noukb_results.txt", header = T) %>% mutate(., overlap = "non-overlapping")

results <- rbind(results_WHR_HNCnoukb,
                 results_BMI_HNCnoukb,
                 results_WC_HNCnoukb)

results$outcome[results$outcome=="OC_noukb"] <- "Oral cavity cancer"
results$outcome[results$outcome=="HPC_noukb"] <- "Hypopharyngeal cancer"
results$outcome[results$outcome=="OPC_POS_noukb"] <- "HPV positive oropharyngeal cancer"
results$outcome[results$outcome=="OPC_NEG_noukb"] <- "HPV negative oropharyngeal cancer"
results$outcome[results$outcome=="LA_noukb"] <- "Laryngeal cancer"
results$outcome[results$outcome=="HNC_noukb"] <- "Head and neck cancer"

#---------------------------------------------------------------------#
#                 Format meta-analysis results                        #----
#---------------------------------------------------------------------#

cancers <- c('Head and neck cancer', 'Oral cavity cancer', 'Laryngeal cancer', 
             'Hypopharyngeal cancer', 'HPV positive oropharyngeal cancer', 'HPV negative oropharyngeal cancer')

as.factor(results$outcome)
results$outcome <- factor(results$outcome, levels = cancers)
results <- arrange(results, outcome)

#reorder methods
as.factor(results$method)
results$method <- factor(results$method, levels = c("Weighted mode", "Weighted median",
                                                    "MR Egger", "Inverse variance weighted"))

results$group[results$outcome!="Head and neck cancer"] <- "Subsites"
results$group[results$outcome=="Head and neck cancer"] <- "Overall"

#---------------------------------------------------------------------#
#                            Plot Function                            #----
#---------------------------------------------------------------------#

myforestplot <- function(df, title, xlab)
{
  options(scipen = 999)
  
  x <- forestplot(
    df = df,
    estimate = b,  
    se = se,
    pvalue = pval,
    name = outcome,
    logodds = TRUE,
    colour = method,
    title = title,
    xlab = xlab,
    #xlim= c(0.88,1.22)
  ) +
    ggforce::facet_col(
      facets = ~group,
      scales = "free_y",
      space = "free"
    )
  colours <- c("#999999", "#E69F00", "#56B4E9", "#009E73")
  x <- x + scale_color_manual(values=colours) 
    print(x)
}

#---------------------------------------------------------------------#
#                            Forest plots                             #----
#---------------------------------------------------------------------#

#main analyses 
BMI_strict_HNC_noukb_plot <- myforestplot(results[results$exposure=="BMI_strict",], "", "Odds ratio (95% CI) per 1-SD higher body mass index")
WHR_strict_HNC_noukb_plot <- myforestplot(results[results$exposure=="WHR_strict",], "", "Odds ratio (95% CI) per 1-SD higher waist-to-hip ratio")
WC_HNC_noukb_plot <- myforestplot(results[results$exposure=="WC",], "", "Odds ratio (95% CI) per 1-SD higher waist circumference")

#---------------------------------------------------------------------#
#                              Save plot                              #----
#---------------------------------------------------------------------#

save_forest_func <- function(file_name, plot_name)
{
  png(file_name, res=330, height=2000, width=3000)
  print(plot_name)
  dev.off()
}

#save 
save_forest_func('results/plots/BMI_strict_HNC_noukb_plot.png', BMI_strict_HNC_noukb_plot)
save_forest_func('results/plots/WHR_strict_HNC_noukb_plot.png', WHR_strict_HNC_noukb_plot)
save_forest_func('results/plots/WC2_HNC_noukb_plot.png', WC_HNC_noukb_plot)

#---------------------------------------------------------------------#
#                              Save csv                              #----
#---------------------------------------------------------------------#

write.csv(results, 'results/MR_results/MR_results_BMI_WHR_WC2.csv', row.names = F)
