#######################################################################
#               HEADSpAcE MENDELIAN RANDOMIZATION ANALYSES            #
#######################################################################
# R version 4.2.3 (2023-03-15)
# Last modified: 23 June 2025

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("YOUR_WD") 

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table")

#---------------------------------------------------------------------#
#                      Power calculations                             #---- 
#---------------------------------------------------------------------#

power_calc_func <- function(cancer_name, Ncases, Ncontrols) {
  
  n <- (Ncases+Ncontrols) 
  ratio <- Ncases/Ncontrols
  
  sig <- 0.05
  
  Betas <- seq(from=0, to=0.5, by=0.0005)
  Powers <- as.data.frame(Betas)
  Powers$ORs <- exp(Betas)
  
  Powers$Var0.5 <- (pnorm(sqrt(n*0.005*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
  Powers$Var1 <- (pnorm(sqrt(n*0.01*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
  Powers$Var2.5 <- (pnorm(sqrt(n*0.025*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
  Powers$Var5 <- (pnorm(sqrt(n*0.05*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
  Powers$Var10 <- (pnorm(sqrt(n*0.1*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
  
  Var4.8 <- (pnorm(sqrt(n*0.048*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
  Var3.1 <- (pnorm(sqrt(n*0.031*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
  Var4.4 <- (pnorm(sqrt(n*0.044*(ratio/(1+ratio))*(1/(1+ratio)))*Powers$Betas-qnorm(1-sig/2)))*100
  
  desired_y_values <- 80
  x_values_approx <- c(approx(Var4.8, Powers$ORs, xout = desired_y_values)$y, approx(Var4.4, Powers$ORs, xout = desired_y_values)$y, approx(Var3.1, Powers$ORs, xout = desired_y_values)$y)
  
  points_data <- data.frame(
    x = x_values_approx,
    y = desired_y_values,
    label=c("BMI", "WC", "WHR"))

  PowerPlot <- ggplot(Powers, aes(ORs)) +
    geom_line(aes(y = Var0.5, colour = "9")) +  
    geom_line(aes(y = Var1, colour = "7")) +
    geom_line(aes(y = Var2.5, colour = "5")) + 
    geom_line(aes(y = Var5, colour = "3")) +
    geom_line(aes(y = Var10, colour = "1")) +
    xlab("Odds ratio per unit increase in risk factor") +
    scale_y_continuous("Power (%)", limits=c(0, 100), breaks=c(20,40,60,80,100)) +
    xlim(1, 1.4) + # Set the x-axis from 1 to 1.4, but make sure data points don't fall out of this range, or else data won't show
    theme(axis.title.x = element_text(face="bold", size=20), axis.text.x  = element_text(vjust=0.5, size=16)) +
    theme(axis.title.y = element_text(face="bold", size=20), axis.text.y  = element_text(vjust=0.5, size=16)) +
    theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
    scale_colour_discrete("% Variance", labels= c("10.0", "5.0", "2.5", "1.0", "0.5")) +
    geom_hline(yintercept=80) +
    ggtitle(paste0("Power for analysing HEADSpAcE ", cancer_name," (", Ncases, " cases & " , Ncontrols, " controls)")) +
    theme(plot.title=element_text(lineheight=5, size= rel(1.2), face="bold")) +
    
    geom_point(data = points_data, aes(x=x, y=y), color="black", size=3, shape=18) +
    geom_text(data = points_data, aes(x=x, y=y, label=label), nudge_y = 3, size=2.5) +
    
    # vertical line and x-axis label
    geom_segment(data = points_data, aes(x = x, y = y, xend = x, yend = 0), 
                 linetype = "dashed", color = "black") +
    geom_text(data = points_data, aes(x = x, y = 0, label = sprintf("%.2f", x)), 
              vjust = 2, color = "black", size = 2.5) +
    
    theme_classic()
  
  png(paste0("results/power_calculations/headspace_", cancer_name, "_power_2SMR.png"), width = 3000, height = 1500, res=200)
  print(PowerPlot)
  dev.off()
  
}

#run function for all cancer outcomes (excluding UKB from cancer outcome)
power_calc_func("HNC", 12264, 19259)
power_calc_func("HPC", 474, 18178)
power_calc_func("LC", 2490, 18178)
power_calc_func("OC", 3091, 18178)
power_calc_func("HPVposOPC", 1980, 18166)
power_calc_func("HPVnegOPC", 948, 18166)

#run function for Gormley et al.
power_calc_func("OC and OPC", 6034, 6585) #I edited the file name and title manually but didn't save


