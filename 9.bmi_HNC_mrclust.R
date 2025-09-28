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
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggplot2", 
               "ggforce", "data.table", "ggforestplot","gtools", "LDlinkR", "simex", "ieugwasr", "openxlsx", "meta", "metafor")
#pacman::p_load_gh("MRCIEU/MRInstruments", "MRCIEU/TwoSampleMR") #use this if you haven't installed these packages before
pacman::p_load_gh("cnfoley/mrclust")

#alternative way of installing and loading MRCIEU packages from GitHub
# install.packages("devtools")
# devtools::install_github("MRCIEU/TwoSampleMR")
# devtools::install_github("MRCIEU/MRInstruments", force = TRUE)
# library(TwoSampleMR)
# library(MRInstruments)

#---------------------------------------------------------------------#
#                                Read MR data                         #----
#---------------------------------------------------------------------#
#BMI Pulit (strict pval)-HNC data
dat <- fread("results/harmonised_data_files/harmonised_BMI_strict_HNC_noukb.csv")

#---------------------------------------------------------------------#
#                               MR analysis                           #----
#---------------------------------------------------------------------#

results <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", 
                                   "mr_weighted_mode",  "mr_egger_regression"))
results

#single snp analysis
res_single <- mr_singlesnp(dat)

#---------------------------------------------------------------------#
#                                  MR-Clust                           #----
#---------------------------------------------------------------------#

res_single <- res_single %>%
  select("SNP", "exposure", "outcome",  "b", "se") %>% 
  filter(grepl("rs", SNP))

dat_sub <- dat %>% 
  select("SNP", "exposure", "outcome",
         "beta.exposure", "se.exposure",
         "beta.outcome", "se.outcome") %>% tidyr::drop_na() 

cluster_data <- inner_join(dat_sub, res_single)

mr_clust_func <- function(outcome_name) {
  
  cluster_data2 <- cluster_data %>% filter(., outcome==outcome_name)
  
  set.seed(123123) #or else results will change every time this is run
  
  cluster_results <- mr_clust_em(theta = cluster_data2$b, 
                                 theta_se = cluster_data2$se,
                                 bx = cluster_data2$beta.exposure, 
                                 by = cluster_data2$beta.outcome,
                                 bxse = cluster_data2$se.exposure,
                                 byse = cluster_data2$se.outcome,
                                 obs_names = cluster_data2$SNP)
  cluster_results$results$best
  clusters <- unique(cluster_results$results$best$cluster_class) 
  print(clusters)
  
  #mean for each cluster
  cluster_means <- unique(cluster_results$results$best$cluster_mean)
  cbind(clusters, cluster_means)
  
  print(cluster_results)
  
}


cluster_res_HNC <- mr_clust_func("HNC_noukb")


#---------------------------------------------------------------------#
#                              Save results                           #----
#---------------------------------------------------------------------#

#save results (best only, as tables)
# write.xlsx(x = cluster_res_HNC$results$best,file = "results/MR_results/mr-clust/mr_clust_BMI_HNC.xlsx")



#save results for the future so we don't have to run them again (you need the large lists for each cancer to run the plots and other analyses)
#save.image(file = "results/MR_results/mr-clust/mr_clust_BMI_HNC.RData")
load("results/MR_results/mr-clust/mr_clust_BMI_HNC.RData")

#create results table
dat_table <- fread("results/harmonised_data_files/harmonised_BMI_strict_HNC_noukb.csv") %>% 
  select("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure") %>% 
  rename(., "observation"= "SNP", "chr"="chr.exposure", "pos"="pos.exposure", "effect_allele"="effect_allele.exposure", "other_allele"="other_allele.exposure")

cluster_res_all <- cluster_res_HNC$results$best %>% mutate(outcome = "HNC") %>% 
  select(observation, outcome, cluster_class, cluster_mean, probability, theta, theta_se) 

col_names_table <- c("outcome", "cluster class","cluster mean", "probability", "beta", "se" )


final_table <- inner_join(dat_table, cluster_res_all)

names(final_table) <- c("SNP", "chromosome", "position", "effect allele", 
                            "other allele", col_names_table)

#write.xlsx(x = final_table,file = "results/MR_results/mr-clust/mr_clust_final_table.xlsx")

#---------------------------------------------------------------------#
#                     Create plots for MR-Clust                       #----
#---------------------------------------------------------------------#

#relevel factors for plot if needed (shouldn't be necessary, especially after filtering the results based on P and nsnps)
#factor(cluster_res_HNC$plots$two_stage$data$clusters, levels = c("3", "2", "1"))
#cluster_res_HNC$plots$two_stage$data <- arrange(cluster_res_HNC$plots$two_stage$data, by = clusters)

#plot results
plot_cluster_func <- function(exp, out, res_object) {
  
  #plot
  clust_plot_best = res_object$plots$two_stage  +
    ggplot2::xlab(paste0("Genetic association with ", exp)) +
    ggplot2::ylab(paste0("Genetic association with ", out)) +
    ggplot2::ggtitle("") # +
    # ggplot2::scale_color_manual(values = c("Null" = "#9944AA22", "Junk" = "#000000", "1" = "#2277DD",
    #                                        "2" = "#E69FFF", "3" = "#D55E00", "4" = "#F0E442",
    #                                        "5" = "#009E73", "6" = "#56B4E9", "7" = "#0072B2")) #no longer necessary since these plots are not the final ones
  
  
  clust_plot_best
  
}

cluster_plot_HNC <- plot_cluster_func("BMI","HNC", cluster_res_HNC)


#keep only best SNPs per cluster (prob > 0.8 and min 4 SNPs per cluster) and rerun plot
plot80_func <- function(exp, out, clust_results) {

res80 = mrclust::pr_clust(dta = clust_results$results$best, prob = 0.8, min_obs =  4)

keep80 = which(clust_results$plots$two_stage$data$observation %in% res80$observation)
bx80   = clust_results$plots$two_stage$data$bx[keep80]
bxse80 = clust_results$plots$two_stage$data$bxse[keep80]
by80   = clust_results$plots$two_stage$data$by[keep80]
byse80 = clust_results$plots$two_stage$data$byse[keep80]
snp_names80 = clust_results$plots$two_stage$data$observation[keep80]

plot.80 = two_stage_plot(res = res80, bx = bx80, by = by80, bxse = bxse80,
                               byse = byse80, obs_names = snp_names80) + 
  ggplot2::xlim(0, max(abs(bx80) + 2*bxse80)) + 
  ggplot2::xlab(paste0("Genetic association with ", exp)) +
  ggplot2::ylab(paste0("Genetic association with ", out)) +
  ggplot2::ggtitle("")

plot.80

}

cluster_plot_HNC80 <- plot80_func("BMI","HNC", cluster_res_HNC)

#side by side plots
combine2_func <- function(A, B){
  plot.new()
  par(mar=c(1,1,1,1), mgp=c(3,1,0))
  x <- ggarrange(A, B, labels=c('A', 'B'),
                 ncol = 2, nrow = 1, common.legend = F, hjust = -3, legend = "right")
  print(x)
}

HNC_plot <- combine2_func(cluster_plot_HNC, cluster_plot_HNC80)

#---------------------------------------------------------------------#
#                                 Save plots                          #----
#---------------------------------------------------------------------#

#save plot
save_func <- function(file_name, plot_name, width)
{
  png(file_name, res=330, height=2000, width=width)
  print(plot_name)
  dev.off()
}


save_func('results/plots/mr-clust_plots/MRclust_plot_BMI_HNC.png', cluster_plot_HNC, 3000)


save_func('results/plots/mr-clust_plots/MRclust_plot_BMI_HNC80.png', cluster_plot_HNC80, 3000)


save_func('results/plots/mr-clust_plots/MRclust_plot_BMI_HNC_combined.png', HNC_plot, 4000)



