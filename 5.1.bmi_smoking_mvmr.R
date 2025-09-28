#######################################################################
#               HEADSpAcE MENDELIAN RANDOMIZATION ANALYSES            #
#######################################################################
# R version 4.4.0 (2024-04-24)
# Last modified: 4 September 2024

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
pacman::p_load_gh("Spiller/MVMR")

#alternative way of installing and loading MRCIEU packages from GitHub
# install.packages("devtools")
# devtools::install_github("MRCIEU/TwoSampleMR")
# devtools::install_github("MRCIEU/MRInstruments", force = TRUE)
# library(TwoSampleMR)
# library(MRInstruments)

#---------------------------------------------------------------------#
#                           Read exposure data                        #----
#---------------------------------------------------------------------#
# #BMI Pulit
# bmi_dat <- fread("results/clumped_exposure_files/clumped_BMI.csv")
# bmi_dat2 <- bmi_dat %>% filter(., pval.exposure <5e-9) #458 snps

# 
# # smoking initiation (si) - excluding UKB data (not the one we need to use)
# si_dat <- fread("results/clumped_exposure_files/clumped_si.csv")
# si_dat2 <- fread("results/clumped_exposure_files/clumped_si_incUKB.csv")
# 
# # comprehensive smoking index
# csi_dat <- fread("results/clumped_exposure_files/clumped_csi.csv")

#---------------------------------------------------------------------#
#                   rsIDs for MVMR analysis                           #----
#---------------------------------------------------------------------#
#read exposure

rsid_func <- function(other_var, full_stats_other_var, snp = "RSID", beta = "BETA", se = "SE", eaf = "AF_1000G", effect_allele = "EFFECT_ALLELE", 
                      other_allele = "OTHER_ALLELE", pval = "P", chr = "CHR", pos = "POS") {
  
  #BMI Pulit
  bmi_dat <- fread("results/clumped_exposure_files/clumped_BMI.csv")
  bmi_dat2 <- bmi_dat %>% filter(., pval.exposure <5e-9) #458 snps

  #other variable (smoking)
  other_exp_dat <- fread(paste0("results/clumped_exposure_files/clumped_", other_var,".csv"))
  other_exp_stats <- fread(full_stats_other_var) %>% data.frame()
  other_exp_stats$pheno <- other_var
  
  #bind datasets
  x <- smartbind(bmi_dat2, other_exp_dat)
  
  #extract list of SNPs that was GWAS-significant in both datasets
  all_snps <-   format_data(snps = x$SNP,
                            other_exp_stats,
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
                            pos_col = pos
  )
  
  
  #remove SNPs in LD across both exposures.
  clumped <- clump_data(all_snps, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1)
  
  #save
  write.csv(clumped, paste0('results/clumped_exposure_files/clumped_mvmr_bmi_', other_var,'.csv'), row.names = F)
  
  print(clumped)
  
}

si_snps <- rsid_func("si_incUKB", "data/GSCAN_2/incUKB/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt.gz") %>% select(-id.exposure)
csi_snps <- rsid_func("csi", "data/comprehensive_smoking_index/csi_completesumstats.txt", snp = "SNP", pos = "BP", eaf = "EAF") %>% select(-id.exposure)


#si_snps <- fread('results/clumped_exposure_files/clumped_mvmr_bmi_si_incUKB.csv') %>% select(-id.exposure)
#csi_snps <- fread('results/clumped_exposure_files/clumped_mvmr_bmi_csi.csv') %>% select(-id.exposure)

#---------------------------------------------------------------------#
#             Read full summary stats for outcomes and bmi            #----
#---------------------------------------------------------------------#

#read stats
bmi_pulit_stats <- fread("data/GIANT_Pulit2018/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")
bmi_pulit_stats$SNP <- gsub("\\:.*","",bmi_pulit_stats$SNP) 

#outcomes
# noukb_HNC <- fread("data/headspace_hnc_gwas_summary_stats_09-2023/noukb_EUR70.nohetrogen_results_header_rsid.txt")
# noukb_HPC <- fread("data/headspace_hnc_gwas_summary_stats_09-2023/noukb_EUR70_HPC.nohetrogen_results_header_rsid.txt")
# noukb_LA <- fread("data/headspace_hnc_gwas_summary_stats_09-2023/noukb_EUR70_LA.nohetrogen_results_header_rsid.txt")
# noukb_OC <- fread("data/headspace_hnc_gwas_summary_stats_09-2023/noukb_EUR70_OC.nohetrogen_results_header_rsid.txt")
# noukb_OPC_NEG <- fread("data/headspace_hnc_gwas_summary_stats_09-2023/noukb_EUR70_OPC_NEG.nohetrogen_results_header_rsid.txt") 
# noukb_OPC_POS <- fread("data/headspace_hnc_gwas_summary_stats_09-2023/noukb_EUR70_OPC_POS.nohetrogen_results_header_rsid.txt")



#---------------------------------------------------------------------#
#                   MVMR exposures and outcomes                       #----
#---------------------------------------------------------------------#

mvmr_prep_func <- function(snps_data, outcome_mrbase, outcome_name) {
  #format other exposure
  other_var_MVMR <- snps_data %>% data.frame()
  # other_var_MVMR <- format_data(other_var_MVMR, type = "exposure", 
  #                        phenotype_col = "outcome",
  #                        snp_col = "SNP",
  #                        beta_col = "beta.outcome", 
  #                        se_col = "se.outcome", 
  #                        eaf_col = "eaf.outcome", 
  #                        effect_allele_col = "effect_allele.outcome", 
  #                        other_allele_col = "other_allele.outcome", 
  #                        pval_col = "pval.outcome", 
  #                        samplesize_col = "samplesize.outcome", 
  #                        chr_col = "chr.outcome", 
  #                        pos_col = "pos.outcome")  
  
  #format BMI exposure
  bmi_stats <- bmi_pulit_stats %>% data.frame()
  bmi_stats$pheno <- "BMI"
  bmi_MVMR <- format_data(snps = snps_data$SNP,
                          bmi_stats,
                          type = "exposure",
                          phenotype_col = "pheno",
                          snp_col = "SNP",
                          beta_col = "BETA",
                          se_col = "SE",
                          eaf_col = "Freq_Tested_Allele",
                          effect_allele_col = "Tested_Allele",
                          other_allele_col = "Other_Allele",
                          pval_col = "P",
                          samplesize_col = "N",
                          chr_col = "CHR",
                          pos_col = "POS"
  )


  #HEADSpAcE outcome
  out_dat <- extract_outcome_data(
    snps = snps_data$SNP,
    outcomes = outcome_mrbase,
    opengwas_jwt = ieugwasr::get_opengwas_jwt()
  )
  out_dat$outcome <- outcome_name


  
  #harmonise
  MVMR_data1 <- harmonise_data(
    exposure_dat = bmi_MVMR,
    outcome_dat = out_dat
  )
  MVMR_data2 <- harmonise_data(
    exposure_dat = bmi_MVMR,
    outcome_dat = other_var_MVMR
  )


  
  data_func <- function(MVMR_DATA) {
    x <- subset(MVMR_DATA, select = c("SNP", "outcome", "exposure","beta.exposure","beta.outcome","se.exposure",
                                      "se.outcome", "pval.exposure", "pval.outcome", "effect_allele.exposure", "effect_allele.outcome", "eaf.exposure", "eaf.outcome"))
  }

  MVMR_data <- merge(data_func(MVMR_data1), data_func(MVMR_data2), by = c("SNP", "exposure", "beta.exposure", "se.exposure", "pval.exposure", "effect_allele.exposure", "eaf.exposure"))

  MVMR_data <- MVMR_data %>% rename(., 
                                    "outcome" = "outcome.x",
                                    "beta.outcome" = "beta.outcome.x",
                                    "se.outcome" = "se.outcome.x",
                                    "pval.outcome" = "pval.outcome.x",
                                    "effect_allele.outcome" = "effect_allele.outcome.x",
                                    "eaf.outcome" = "eaf.outcome.x",
                                    "exposure.y" = "outcome.y",
                                    "beta.exposure.y" = "beta.outcome.y",
                                    "se.exposure.y" = "se.outcome.y",
                                    "pval.exposure.y" = "pval.outcome.y",
                                    "effect_allele.exposure.y" = "effect_allele.outcome.y",
                                    "eaf.exposure.y" = "eaf.outcome.y",
                                    "exposure.x" = "exposure",
                                    "beta.exposure.x" = "beta.exposure",
                                    "se.exposure.x" = "se.exposure",
                                    "pval.exposure.x" = "pval.exposure",
                                    "effect_allele.exposure.x" = "effect_allele.exposure",
                                    "eaf.exposure.x" = "eaf.exposure")
  
}


#si
bmi_si_HNC_data <- mvmr_prep_func(si_snps, 'ieu-b-5123',"HNC_noukb")
bmi_si_HPC_data <- mvmr_prep_func(si_snps, 'ieu-b-5124',"HPC_noukb")
bmi_si_LA_data <- mvmr_prep_func(si_snps, 'ieu-b-5125',"LA_noukb")
bmi_si_OC_data <- mvmr_prep_func(si_snps, 'ieu-b-5126',"OC_noukb")
bmi_si_OPC_NEG_data <- mvmr_prep_func(si_snps, 'ieu-b-5127',"OPC_NEG_noukb")
bmi_si_OPC_POS_data <- mvmr_prep_func(si_snps, 'ieu-b-5128',"OPC_POS_noukb")

#csi
bmi_csi_HNC_data <- mvmr_prep_func(csi_snps, 'ieu-b-5123',"HNC_noukb")
bmi_csi_HPC_data <- mvmr_prep_func(csi_snps, 'ieu-b-5124',"HPC_noukb")
bmi_csi_LA_data <- mvmr_prep_func(csi_snps, 'ieu-b-5125',"LA_noukb")
bmi_csi_OC_data <- mvmr_prep_func(csi_snps, 'ieu-b-5126',"OC_noukb")
bmi_csi_OPC_NEG_data <- mvmr_prep_func(csi_snps, 'ieu-b-5127',"OPC_NEG_noukb")
bmi_csi_OPC_POS_data <- mvmr_prep_func(csi_snps, 'ieu-b-5128',"OPC_POS_noukb")


#---------------------------------------------------------------------#
#                    save harmonised data                             #----
#---------------------------------------------------------------------#
#data was harmonised on the same strand (using same exposure for both harmonisation processes)

#save
save_harmonised_func <- function(harmonised_data, other_var, outcome) {
  write.csv(harmonised_data, paste0('results/harmonised_data_files/mvmr/harmonised_mvmr_bmi_', other_var, '_', outcome, '.csv'), row.names = F)
}

save_harmonised_func(bmi_si_HNC_data, "SI_incUKB", "HNC_noukb")
save_harmonised_func(bmi_si_HPC_data, "SI_incUKB", "HPC_noukb")
save_harmonised_func(bmi_si_LA_data, "SI_incUKB", "LA_noukb")
save_harmonised_func(bmi_si_OC_data, "SI_incUKB", "OC_noukb")
save_harmonised_func(bmi_si_OPC_NEG_data, "SI_incUKB", "OPC_NEG_noukb")
save_harmonised_func(bmi_si_OPC_POS_data, "SI_incUKB", "OPC_POS_noukb")

save_harmonised_func(bmi_csi_HNC_data, "CSI", "HNC_noukb")
save_harmonised_func(bmi_csi_HPC_data, "CSI", "HPC_noukb")
save_harmonised_func(bmi_csi_LA_data, "CSI", "LA_noukb")
save_harmonised_func(bmi_csi_OC_data, "CSI", "OC_noukb")
save_harmonised_func(bmi_csi_OPC_NEG_data, "CSI", "OPC_NEG_noukb")
save_harmonised_func(bmi_csi_OPC_POS_data, "CSI", "OPC_POS_noukb")


#---------------------------------------------------------------------#
#                           MVMR analysis                             #----
#---------------------------------------------------------------------#

#MVMR function
mvmr_func <- function(data, outcome) {
  x <- format_mvmr(BXGs = data[,c("beta.exposure.x", "beta.exposure.y")], BYG = data[,"beta.outcome"], seBXGs = data[,c("se.exposure.x", "se.exposure.y")], seBYG = data[,"se.outcome"], RSID = data[,"SNP"])
  a <- as.data.frame(t(strength_mvmr(x, gencov = 0)), row.names = c(paste0(data$exposure.x[1], "_MVMR", sep = ""), paste0(data$exposure.y[1], "_MVMR", sep = "")))
  y <- as.data.frame(ivw_mvmr(x, gencov = 0), row.names = c(paste0(data$exposure.x[1], "_MVMR", sep = ""), paste0(data$exposure.y[1], "_MVMR", sep = "")))
  y$outcome <- outcome
  y$exposure.1 <- data$exposure.x[1]
  y$exposure.2 <- data$exposure.y[1]
  m <- cbind(y, a)
  return(m)
}


#Run MVMR
#si
mvmr_bmi_si_HNC_result <- mvmr_func(bmi_si_HNC_data, "HNC")
mvmr_bmi_si_HPC_result <- mvmr_func(bmi_si_HPC_data, "HPC")
mvmr_bmi_si_LA_result <- mvmr_func(bmi_si_LA_data, "LA")
mvmr_bmi_si_OC_result <- mvmr_func(bmi_si_OC_data, "OC")
mvmr_bmi_si_OPC_NEG_result <- mvmr_func(bmi_si_OPC_NEG_data, "OPC_NEG")
mvmr_bmi_si_OPC_POS_result <- mvmr_func(bmi_si_OPC_POS_data, "OPC_POS")


results_bmi_si <- rbind(mvmr_bmi_si_HNC_result, 
                     mvmr_bmi_si_HPC_result, 
                     mvmr_bmi_si_LA_result, 
                     mvmr_bmi_si_OC_result, 
                     mvmr_bmi_si_OPC_NEG_result,
                     mvmr_bmi_si_OPC_POS_result)



#csi
mvmr_bmi_csi_HNC_result <- mvmr_func(bmi_csi_HNC_data, "HNC")
mvmr_bmi_csi_HPC_result <- mvmr_func(bmi_csi_HPC_data, "HPC")
mvmr_bmi_csi_LA_result <- mvmr_func(bmi_csi_LA_data, "LA")
mvmr_bmi_csi_OC_result <- mvmr_func(bmi_csi_OC_data, "OC")
mvmr_bmi_csi_OPC_NEG_result <- mvmr_func(bmi_csi_OPC_NEG_data, "OPC_NEG")
mvmr_bmi_csi_OPC_POS_result <- mvmr_func(bmi_csi_OPC_POS_data, "OPC_POS")


results_bmi_csi <- rbind(mvmr_bmi_csi_HNC_result, 
                        mvmr_bmi_csi_HPC_result, 
                        mvmr_bmi_csi_LA_result, 
                        mvmr_bmi_csi_OC_result, 
                        mvmr_bmi_csi_OPC_NEG_result,
                        mvmr_bmi_csi_OPC_POS_result) #remember to rescale the MVMR results before creating plots if CSI results adjusted for BMI are needed

#csi per SD increase if this hasn't been done above on "dat"
results_bmi_csi[str_detect(row.names(results_bmi_csi), "csi"),]$Estimate <- results_bmi_csi[str_detect(row.names(results_bmi_csi), "csi"),]$Estimate*0.6940093
results_bmi_csi[str_detect(row.names(results_bmi_csi), "csi"),]$Std..Error <- results_bmi_csi[str_detect(row.names(results_bmi_csi), "csi"),]$Std..Error*0.6940093

#save results
write.xlsx(results_bmi_si, "results/MVMR_results/MVMR_BMI_SI_incUKB.xlsx", row.names= T)
write.xlsx(results_bmi_csi, "results/MVMR_results/MVMR_BMI_CSI.xlsx", row.names= T)

#results_bmi_csi <- read.xlsx("results/MVMR_results/MVMR_BMI_CSI.xlsx", rowNames = T)
#results_bmi_si <- read.xlsx("results/MVMR_results/MVMR_BMI_SI_incUKB.xlsx", rowNames = T)

#################################################################################


