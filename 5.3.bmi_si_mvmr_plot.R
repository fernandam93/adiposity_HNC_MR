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

if (!require("pacman")) install.packages("pacman")
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot", "gtools", "openxlsx", "remotes", "purrr")

#---------------------------------------------------------------------#
#                            MVMR data prep                           #----
#---------------------------------------------------------------------#
##LOAD MVMR DATA
MVMR_BMI_SI <- read.xlsx("results/MVMR_results/MVMR_BMI_SI_incUKB.xlsx", rowNames = T)

names(MVMR_BMI_SI)[names(MVMR_BMI_SI)=="Estimate"] <- "b"
names(MVMR_BMI_SI)[names(MVMR_BMI_SI)=="Std..Error"] <- "se"
names(MVMR_BMI_SI)[names(MVMR_BMI_SI)=="Pr(>|t|)"] <- "pval"
names(MVMR_BMI_SI)[names(MVMR_BMI_SI)=="exposure.1"] <- "exposure"
MVMR_BMI_SI <- MVMR_BMI_SI[grepl("BMI", rownames(MVMR_BMI_SI)), ]

MVMR_BMI_SI$method <- "Multivariable MR"

MVMR_BMI_SI <- subset(MVMR_BMI_SI, select = -c(t.value, `F-statistic`, exposure.2))

#---------------------------------------------------------------------#
#                            Univar data prep                         #----
#---------------------------------------------------------------------#
##LOAD UNIVAR DATA
univar_BMI <- fread("results/MR_results/BMI_strict_hnc_noukb_results.txt")
univar_BMI <- univar_BMI[univar_BMI$exposure=="BMI_strict" & univar_BMI$method=="Inverse variance weighted",]

univar_BMI$method <- "Univariable MR"

univar_BMI <- subset(univar_BMI, select = -c(or, or_lci95, or_uci95, nsnp))

univar_BMI$exposure <- "BMI"

#---------------------------------------------------------------------#
#                       Merge univar and MVMR datasets                #----
#---------------------------------------------------------------------#

results_BMI <- rbind(univar_BMI, MVMR_BMI_SI)

results_BMI$outcome[results_BMI$outcome=="OC_noukb"] <- "Oral cavity cancer"
results_BMI$outcome[results_BMI$outcome=="HPC_noukb"] <- "Hypopharyngeal cancer"
results_BMI$outcome[results_BMI$outcome=="OPC_POS_noukb"] <- "HPV positive oropharyngeal cancer"
results_BMI$outcome[results_BMI$outcome=="OPC_NEG_noukb"] <- "HPV negative oropharyngeal cancer"
results_BMI$outcome[results_BMI$outcome=="LA_noukb"] <- "Laryngeal cancer"
results_BMI$outcome[results_BMI$outcome=="HNC_noukb"] <- "Head and neck cancer"


#---------------------------------------------------------------------#
#                           Plot preparation                          #----
#---------------------------------------------------------------------#

cancers <- c('Head and neck cancer', 'Oral cavity cancer', 'Laryngeal cancer', 
             'Hypopharyngeal cancer', 'HPV positive oropharyngeal cancer', 'HPV negative oropharyngeal cancer')

as.factor(results_BMI$outcome)
results_BMI$outcome <- factor(results_BMI$outcome, levels = cancers)
results_BMI <- arrange(results_BMI, outcome)

as.factor(results_BMI$method)
results_BMI$method <- factor(results_BMI$method, levels = c("Multivariable MR", "Univariable MR"))

results_BMI$group[results_BMI$outcome!="Head and neck cancer"] <- "Subsites"
results_BMI$group[results_BMI$outcome=="Head and neck cancer"] <- "Overall"



#---------------------------------------------------------------------#
#                        Forest plot function                         #----
#---------------------------------------------------------------------#
#BINARY OUTCOMES
myforestplot_mvmr <- function(df, title, xlab)
{
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
    #xlim= c(0,110)
  ) +
    ggforce::facet_col(
      facets = ~group,
      scales = "free_y",
      space = "free"
    )
  colours_BP <- c("#3899FF", "black")
  x <- x + scale_color_manual(values=colours_BP)
  x <- x + scale_x_continuous(breaks = c(0.3, 0.5, 0.7, 1, 1.3, 1.5, 1.7))
  print(x)
}


#---------------------------------------------------------------------#
#                                Forest plot                          #----
#---------------------------------------------------------------------#

p_mvmr <- myforestplot_mvmr(results_BMI, "", "OR (95% CI) per 1-SD higher body mass index")

#---------------------------------------------------------------------#
#                                 Save plot                           #----
#---------------------------------------------------------------------#

save_forest_func <- function(file_name, plot_name)
{
  png(file_name, res=330, height=2000, width=3000)
  print(plot_name)
  dev.off()
}

save_forest_func("results/plots/mvmr_plots/mvmr_bmi_si_incUKB.png", p_mvmr)

