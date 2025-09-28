######################################################################
##                MR-CAUSE - BMI and Head and Neck Cancer              ##
######################################################################

# R version 4.4.0 (2024-04-24)
# Last modified: 24 July 2024

# Very slow, so it's better to use Linux

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory - folder in my computer
setwd("YOUR_WD") 

if (!require("pacman")) install.packages("pacman")
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "readr", "ieugwasr", "data.table")
pacman::p_load_gh("jean997/cause@v1.2.0", "explodecomputer/genetics.binaRies")
#devtools::install_github("jean997/cause@v1.2.0")

#---------------------------------------------------------------------#
#                            Read exposure                             #----
#---------------------------------------------------------------------#
#read exposure
bmi_stats <- fread("bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")
bmi_stats$SNP <- gsub("\\:.*","",bmi_stats$SNP)

#---------------------------------------------------------------------#
#                            Outcomes                                  #----
#---------------------------------------------------------------------#
#IEU catalog (MR-Base) - HEADSpAcE data (exceeds max)
HNC <- fread("noukb_EUR70.nohetrogen_results_header_rsid.txt.gz", header = T)

#---------------------------------------------------------------------#
#                            Merge GWAS data                          #----
#---------------------------------------------------------------------#

X <- gwas_merge(bmi_stats, HNC, 
                snp_name_cols = c("SNP", "RSID"), 
                beta_hat_cols = c("BETA", "BETA"), 
                se_cols = c("SE", "SE"), 
                A1_cols = c("Tested_Allele", "A1"), 
                A2_cols = c("Other_Allele", "AX"))

#---------------------------------------------------------------------#
#                    Calculate nuisance parameters                    #----
#---------------------------------------------------------------------#
# only > 100,000 variants
set.seed(100)
varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X, varlist)
head(params$mix_grid)

#---------------------------------------------------------------------#
#                                Clump data                           #----
#---------------------------------------------------------------------#
X$p_value <- 2*pnorm(abs(X$beta_hat_1/X$seb1), lower.tail=FALSE)
X_clump <- X %>% rename(rsid = snp,
                        pval = p_value) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = 0.001,
                     clump_p = 1e-03,
                     # plink_bin = genetics.binaRies::get_plink_binary(),
                     #bfile = "~/EUR"
  )
keep_snps <- X_clump$rsid
#---------------------------------------------------------------------#
#                    MR-CAUSE analysis                                #----
#---------------------------------------------------------------------#
# X unclumped data and variants clumped data
res <- cause(X=X, variants = keep_snps, param_ests = params)
plot(res$sharing)
plot(res$causal)
summary(res, ci_size=0.95)
plot(res)
plot(res, type="data")

png('results/CAUSE_BMI_HNC1.png', res=300, height=2000, width=3500)
plot(res)
dev.off()

png('results/CAUSE_BMI_HNC2.png', res=300, height=2000, width=3500)
plot(res, type="data")
dev.off()
