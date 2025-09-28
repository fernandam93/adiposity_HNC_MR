#######################################################################
#               HEADSpAcE MENDELIAN RANDOMIZATION ANALYSES            #
#######################################################################
# R version 4.4.0 (2024-04-24)
# Last modified: 24 June 2024

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


#---------------------------------------------------------------------#
#                           Read harmonised data                        #----
#---------------------------------------------------------------------#

BMI_HNCnoukb <- read.csv("results/harmonised_data_files/harmonised_BMI_strict_HNC_noukb.csv", header = T) 
BMI_HPCnoukb <- read.csv("results/harmonised_data_files/harmonised_BMI_strict_HPC_noukb.csv", header = T) 
BMI_OCnoukb <- read.csv("results/harmonised_data_files/harmonised_BMI_strict_OC_noukb.csv", header = T) 
BMI_LAnoukb <- read.csv("results/harmonised_data_files/harmonised_BMI_strict_LA_noukb.csv", header = T) 
BMI_OPC_NEGnoukb <- read.csv("results/harmonised_data_files/harmonised_BMI_strict_OPC_NEG_noukb.csv", header = T) 
BMI_OPC_POSnoukb <- read.csv("results/harmonised_data_files/harmonised_BMI_strict_OPC_POS_noukb.csv", header = T) 

#---------------------------------------------------------------------#
#                  Steiger filtering confounder prep                  #----
#---------------------------------------------------------------------#

BMI <-  BMI_HNCnoukb %>% select(SNP, exposure, id.exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure, samplesize.exposure, chr.exposure, pos.exposure)

harmonise_exp_conf_func <- function(other_exp_stats, other_var, snp = "RSID", beta = "BETA", se = "SE", eaf = "AF_1000G", effect_allele = "EFFECT_ALLELE", 
                                    other_allele = "OTHER_ALLELE", pval = "P", chr = "CHR", pos = "POS") {
  
  #we need to create harmonised datasets for BMI-smoking
  #other variable (smoking)
  smoking_data <- fread(other_exp_stats) %>% data.frame()
  smoking_data$pheno <- other_var
  
  smoking <- format_data(snps = BMI$SNP,
                         smoking_data,
                         type = "outcome",
                         phenotype_col = "pheno",
                         snp_col = snp,
                         beta_col = beta,
                         se_col = se,
                         eaf_col = eaf,
                         effect_allele_col = effect_allele,
                         other_allele_col = other_allele,
                         pval_col = pval,
                         samplesize_col = "N",
                         chr_col = chr,
                         pos_col = pos)
  
  dat <- harmonise_data(
    exposure_dat = BMI,
    outcome_dat = smoking)
  
}

bmi_si <- harmonise_exp_conf_func("data/GSCAN_2/incUKB/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt.gz", "si_incUKB")
bmi_csi <- harmonise_exp_conf_func("data/comprehensive_smoking_index/csi_completesumstats.txt", "csi", snp = "SNP", pos = "BP", eaf = "EAF")

bmi_csi$samplesize.outcome <- 462690

#steiger filtering
dat2_si <- steiger_filtering(bmi_si)
dat2_csi <- steiger_filtering(bmi_csi)
#list SNPs that do not pass steiger test
dat3_si <- dat2_si[which(dat2_si$steiger_dir=="FALSE"),]
dat3_csi <- dat2_csi[which(dat2_csi$steiger_dir=="FALSE"),]
#only select SNPs that passed steiger filtering for both smoking traits
dat4a <- BMI_HNCnoukb[which(!(BMI_HNCnoukb$SNP %in% dat3_si$SNP) & !(BMI_HNCnoukb$SNP %in% dat3_csi$SNP)),]
dat4b <- BMI_HPCnoukb[which(!(BMI_HPCnoukb$SNP %in% dat3_si$SNP) & !(BMI_HPCnoukb$SNP %in% dat3_csi$SNP)),]
dat4c <- BMI_OCnoukb[which(!(BMI_OCnoukb$SNP %in% dat3_si$SNP) & !(BMI_OCnoukb$SNP %in% dat3_csi$SNP)),]
dat4d <- BMI_LAnoukb[which(!(BMI_LAnoukb$SNP %in% dat3_si$SNP) & !(BMI_LAnoukb$SNP %in% dat3_csi$SNP)),]
dat4e <- BMI_OPC_NEGnoukb[which(!(BMI_OPC_NEGnoukb$SNP %in% dat3_si$SNP) & !(BMI_OPC_NEGnoukb$SNP %in% dat3_csi$SNP)),]
dat4f <- BMI_OPC_POSnoukb[which(!(BMI_OPC_POSnoukb$SNP %in% dat3_si$SNP) & !(BMI_OPC_POSnoukb$SNP %in% dat3_csi$SNP)),]

dat <- rbind(dat4a, dat4b, dat4c, dat4d, dat4e, dat4f)

#rerun analyses using filtered data
results <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", 
                                   "mr_weighted_median", "mr_weighted_mode"))


write.table(results, "results/MR_results/BMI_strict_hnc_noukb_steiger_results.txt", row.names = F) 
