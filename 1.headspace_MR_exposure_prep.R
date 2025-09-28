#######################################################################
#               HEADSpAcE MENDELIAN RANDOMIZATION ANALYSES            #
#######################################################################
# R version 4.2.3 (2023-03-15)

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("YOUR_WD") 

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", 
               "data.table", "ggforestplot","gtools", "LDlinkR", "simex", "ieugwasr", "R.utils", "openxlsx")
#pacman::p_load_gh("MRCIEU/MRInstruments", "MRCIEU/TwoSampleMR") #use this if you haven't installed these packages before

#alternative way of installing and loading MRCIEU packages from GitHub
# install.packages("devtools")
# devtools::install_github("MRCIEU/TwoSampleMR")
# devtools::install_github("MRCIEU/MRInstruments", force = TRUE)
# library(TwoSampleMR)
# library(MRInstruments)
# install.packages('R.utils')

#---------------------------------------------------------------------#
#                   Prepare exposure data for analysis                #----
#---------------------------------------------------------------------#

#0. pre-process GWAS file if it's too big

#BMI
bmi_stats <- fread("data/GIANT_Pulit2018/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz") %>%
  filter(., P < 5e-08) #I will use a stricter threshold in my analyses, but decided to use the traditional threshold at this stage (will do it in the next script)
bmi_stats$SNP <- gsub("\\:.*","",bmi_stats$SNP)
fwrite(x = bmi_stats, file = "data/GIANT_Pulit2018/bmi.giant-ukbb.smallpval.txt.gz", sep = "\t",quote = F, compress = "gzip") 

#WHR
whr_stats <- fread("data/GIANT_Pulit2018/whr.giant-ukbb.meta-analysis.combined.23May2018.txt.gz") %>%
  filter(., P < 5e-08) #I will use a stricter threshold in my analyses, but decided to use the traditional threshold at this stage (will do it in the next script)
whr_stats$SNP <- gsub("\\:.*","",whr_stats$SNP)
fwrite(x = whr_stats, file = "data/GIANT_Pulit2018/whr.giant-ukbb.smallpval.txt.gz", sep = "\t",quote = F, compress = "gzip")

#waist circumference
wc_stats2 <- extract_instruments(
  c("ukb-b-9405"),
  p1 = 5e-08,
  clump = F,
  opengwas_jwt = ieugwasr::get_opengwas_jwt(),
  force_server = FALSE
)
fwrite(x = wc_stats2, file = "data/UKB/wc.ukb2018elsworth.smallpval.txt.gz", sep = "\t",quote = F, compress = "gzip")

#favourable adiposity
fa_stats <- read.xlsx(xlsxFile = "data/FA_UFA_Martin2021/UFA_FA.xlsx", sheet = "FA", colNames = T) %>% 
  filter(., P < 5e-08)
fwrite(x = fa_stats, file = "data/FA_UFA_Martin2021/FA.smallpval.txt.gz", sep = "\t",quote = F, compress = "gzip")

#unfavourable adiposity
ufa_stats <- read.xlsx(xlsxFile = "data/FA_UFA_Martin2021/UFA_FA.xlsx", sheet = "UFA", colNames = T) %>% 
  filter(., P < 5e-08)
fwrite(x = ufa_stats, file = "data/FA_UFA_Martin2021/UFA.smallpval.txt.gz", sep = "\t",quote = F, compress = "gzip")

#principal components of body shape
pc1_stats <- fread("data/GIANT_bodyshape_Ried2016/GIANT_metal_result_bodyshape_pc1_all_iv_hetero_20111006_adjusted1.txt.gz") %>%
  filter(., p < 5e-08)
pc1_stats <- pc1_stats[!is.na(pc1_stats$CHR_B37)]
fwrite(x = pc1_stats, file = "data/GIANT_bodyshape_Ried2016/GIANT_metal_result_bodyshape_pc1.smallpval.txt.gz", sep = "\t",quote = F, compress = "gzip")

pc2_stats <- fread("data/GIANT_bodyshape_Ried2016/GIANT_metal_result_bodyshape_pc2_all_iv_hetero_20111006_adjusted1.txt.gz") %>%
  filter(., p < 5e-08)
pc2_stats <- pc2_stats[!is.na(pc2_stats$CHR_B37)]
fwrite(x = pc2_stats, file = "data/GIANT_bodyshape_Ried2016/GIANT_metal_result_bodyshape_pc2.smallpval.txt.gz", sep = "\t",quote = F, compress = "gzip")

pc3_stats <- fread("data/GIANT_bodyshape_Ried2016/GIANT_metal_result_bodyshape_pc3_all_iv_hetero_20111006_adjusted1.txt.gz") %>%
  filter(., p < 5e-08)
pc3_stats <- pc3_stats[!is.na(pc3_stats$CHR_B37)]
fwrite(x = pc3_stats, file = "data/GIANT_bodyshape_Ried2016/GIANT_metal_result_bodyshape_pc3.smallpval.txt.gz", sep = "\t",quote = F, compress = "gzip")

pc4_stats <- fread("data/GIANT_bodyshape_Ried2016/GIANT_metal_result_bodyshape_pc4_all_iv_hetero_20111006_adjusted1.txt.gz") %>%
  filter(., p < 5e-08)
pc4_stats <- pc4_stats[!is.na(pc4_stats$CHR_B37)]
fwrite(x = pc4_stats, file = "data/GIANT_bodyshape_Ried2016/GIANT_metal_result_bodyshape_pc4.smallpval.txt.gz", sep = "\t",quote = F, compress = "gzip")

#smoking initiation
si_stats2 <- fread("data/GSCAN_2/incUKB/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt.gz")%>%
  filter(., P < 5e-08) 
fwrite(x = si_stats2, file = "data/GSCAN_2/incUKB/si_incUKB.smallpval.txt.gz", sep = "\t",quote = F, compress = "gzip")


#1. Extract SNP summary stats for adiposity from GWAS and format for MR analysis

exposure_prep_func <- function(exposure_input_file, exposure_name,   
                               snp="SNP", beta="BETA", se="SE", 
                               eaf="FREQ1", effect_allele="A1",other_allele="A0",
                               pval="P", chr="CHR", pos="POS", sep = "\t") {
  
  #read and clump data
  exposure_stats <- read_exposure_data(
    exposure_input_file,
    clump = TRUE,
    sep = sep,
    snp_col = snp,
    beta_col = beta,
    se_col = se,
    eaf_col = eaf,
    effect_allele_col = effect_allele,
    other_allele_col = other_allele,
    pval_col = pval,
    min_pval = 1e-200,
    log_pval = FALSE,
    chr_col = chr,
    pos_col = pos, 
    samplesize_col = "N"
  )
  exposure_stats$exposure <- exposure_name
  
  #save
  write.csv(exposure_stats, paste0('results/clumped_exposure_files/clumped_', exposure_name,'.csv'), row.names = F)
  
  exposure_stats
  
}

#apply function to adiposity traits

bmi_dat <- exposure_prep_func("data/GIANT_Pulit2018/bmi.giant-ukbb.smallpval.txt.gz", 
                              "BMI", eaf = "Freq_Tested_Allele", effect_allele = "Tested_Allele", 
                              other_allele = "Other_Allele")

whr_dat <- exposure_prep_func("data/GIANT_Pulit2018/whr.giant-ukbb.smallpval.txt.gz", 
                              "WHR", eaf = "Freq_Tested_Allele", effect_allele = "Tested_Allele", 
                              other_allele = "Other_Allele")

height_dat <- exposure_prep_func("data/GIANT_height_Yengo2022/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.smallpval.txt.gz", 
                              "Height", eaf = "EFFECT_ALLELE_FREQ", effect_allele = "EFFECT_ALLELE", 
                              other_allele = "OTHER_ALLELE", snp = "RSID")

wc_dat2 <- exposure_prep_func("data/UKB/wc.ukb2018elsworth.smallpval.txt.gz", 
                             "WC2", beta = "beta.exposure", eaf = "eaf.exposure", effect_allele = "effect_allele.exposure", 
                             other_allele = "other_allele.exposure", se = "se.exposure", pval = "pval.exposure",
                             chr = "chr.exposure", pos = "pos.exposure")

fa_dat <- exposure_prep_func("data/FA_UFA_Martin2021/FA.smallpval.txt.gz", 
                                 "FA", eaf = "Freq", effect_allele = "Effect.allele", 
                                 other_allele = "Other.allele", snp = "RSID", chr = "Chromosome", pos = "Position", beta = "Beta")

ufa_dat <- exposure_prep_func("data/FA_UFA_Martin2021/UFA.smallpval.txt.gz", 
                             "UFA", eaf = "Freq", effect_allele = "Effect.allele", 
                             other_allele = "Other.allele", snp = "RSID", chr = "Chromosome", pos = "Position", beta = "Beta")

pc1_dat <- exposure_prep_func("data/GIANT_bodyshape_Ried2016/GIANT_metal_result_bodyshape_pc1.smallpval.txt.gz", 
                              "PC1", eaf = "Freq_Allele1_HapMapCEU", effect_allele = "Allele1", 
                              other_allele = "Allele2", snp = "MarkerName", beta = "beta", pval = "p", chr = "CHR_B37", pos = "POS_B37")

pc2_dat <- exposure_prep_func("data/GIANT_bodyshape_Ried2016/GIANT_metal_result_bodyshape_pc2.smallpval.txt.gz", 
                              "PC2", eaf = "Freq_Allele1_HapMapCEU", effect_allele = "Allele1", 
                              other_allele = "Allele2", snp = "MarkerName", beta = "beta", pval = "p", chr = "CHR_B37", pos = "POS_B37")
pc3_dat <- exposure_prep_func("data/GIANT_bodyshape_Ried2016/GIANT_metal_result_bodyshape_pc3.smallpval.txt.gz", 
                              "PC3", eaf = "Freq_Allele1_HapMapCEU", effect_allele = "Allele1", 
                              other_allele = "Allele2", snp = "MarkerName", beta = "beta", pval = "p", chr = "CHR_B37", pos = "POS_B37")
pc4_dat <- exposure_prep_func("data/GIANT_bodyshape_Ried2016/GIANT_metal_result_bodyshape_pc4.smallpval.txt.gz", 
                              "PC4", eaf = "Freq_Allele1_HapMapCEU", effect_allele = "Allele1", 
                              other_allele = "Allele2", snp = "MarkerName", beta = "beta", pval = "p", chr = "CHR_B37", pos = "POS_B37")


#Smoking traits used in MVMR


csi_dat <- exposure_prep_func("data/comprehensive_smoking_index/csi_exposure.txt", 
                              "csi", beta = "beta", se = "se", eaf = "eaf", effect_allele = "effect_allele", 
                              other_allele = "other_allele", pval = "pval")

si_dat2 <- exposure_prep_func("data/GSCAN_2/incUKB/si_incUKB.smallpval.txt.gz", 
                             "si_incUKB", snp = "RSID", beta = "BETA", se = "SE", eaf = "AF_1000G", effect_allele = "EFFECT_ALLELE", 
                             other_allele = "OTHER_ALLELE", pval = "P")



