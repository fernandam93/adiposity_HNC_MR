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
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot","gtools", "LDlinkR", "simex", "ieugwasr", "openxlsx")
#pacman::p_load_gh("MRCIEU/MRInstruments", "MRCIEU/TwoSampleMR") #use this if you haven't installed these packages before

#---------------------------------------------------------------------#
#                             Read harmonised data                    #----
#---------------------------------------------------------------------#

dat_BMI_strict_HNC <- fread("results/harmonised_data_files/harmonised_BMI_strict_HNC_noukb.csv", header = T) 
dat_WHR_strict_HNC <- fread("results/harmonised_data_files/harmonised_WHR_strict_HNC_noukb.csv", header = T) 
dat_WC_HNC <- fread("results/harmonised_data_files/harmonised_WC_HNC_noukb.csv", header = T) 

dat_adipose_HNC <- fread("results/harmonised_data_files/harmonised_adipose-BMI_HNC_noukb.csv", header = T) 
dat_brain_HNC <- fread("results/harmonised_data_files/harmonised_brain-BMI_HNC_noukb.csv", header = T) 
dat_adult_HNC <- fread("results/harmonised_data_files/harmonised_adult BMI_HNC_noukb.csv", header = T) 
dat_child_HNC <- fread("results/harmonised_data_files/harmonised_child BMI_HNC_noukb.csv", header = T) 
dat_BF_HNC <- fread("results/harmonised_data_files/harmonised_Body fat percentage_HNC_noukb.csv", header = T) 
dat_FA_HNC <- fread("results/harmonised_data_files/harmonised_FA_HNC_noukb.csv", header = T) 
dat_UFA_HNC <- fread("results/harmonised_data_files/harmonised_UFA_HNC_noukb.csv", header = T) 
dat_PC1_HNC <- fread("results/harmonised_data_files/harmonised_PC1_HNC_noukb.csv", header = T) 
dat_PC2_HNC <- fread("results/harmonised_data_files/harmonised_PC2_HNC_noukb.csv", header = T) 
dat_PC3_HNC <- fread("results/harmonised_data_files/harmonised_PC3_HNC_noukb.csv", header = T) 
dat_PC4_HNC <- fread("results/harmonised_data_files/harmonised_PC4_HNC_noukb.csv", header = T) 


dat <- smartbind(dat_BMI_strict_HNC, dat_WHR_strict_HNC, dat_WC_HNC)

dat$samplesize.outcome <- 31523 #total HEADSpAcE sample size (31,523 for sample excluding UKB, 46,468 for the sample including UKB)
dat$samplesize.exposure[dat$exposure=="WC"] <- 462166 
##dat$units.outcome <- "log odds"
##dat$units.exposure <- "SD"

dat <- steiger_filtering(dat) #this is to obtain r2 for the exposure
dat$F <- dat$rsq.exposure * (dat$samplesize.exposure - 2) / (1 - dat$rsq.exposure)


#secondary analyses
dat_secondary <- smartbind(dat_adipose_HNC, dat_brain_HNC, dat_adult_HNC, dat_child_HNC,
                           dat_BF_HNC, dat_FA_HNC, dat_UFA_HNC, dat_PC1_HNC,
                           dat_PC2_HNC, dat_PC3_HNC, dat_PC4_HNC)
dat_secondary$samplesize.outcome <- 31523 #total HEADSpAcE sample size
dat_secondary$samplesize.exposure[dat_secondary$exposure=="FA"] <- 451099 
dat_secondary$samplesize.exposure[dat_secondary$exposure=="UFA"] <- 451099 
dat_secondary$samplesize.exposure[dat_secondary$exposure=="brain-BMI"] <- 681275 
dat_secondary$samplesize.exposure[dat_secondary$exposure=="adipose-BMI"] <- 681275 
dat_secondary$samplesize.exposure[dat_secondary$exposure=="child BMI"] <- 453169 
dat_secondary$samplesize.exposure[dat_secondary$exposure=="adult BMI"] <- 453169 
dat_secondary$samplesize.exposure <- as.numeric(dat_secondary$samplesize.exposure)

dat_secondary <- steiger_filtering(dat_secondary) #this is to obtain r2 for the exposure
dat_secondary$F <- dat_secondary$rsq.exposure * (dat_secondary$samplesize.exposure - 2) / (1 - dat_secondary$rsq.exposure)

dat <- rbind(dat, dat_secondary)

#---------------------------------------------------------------------#
#                       r-squared and F-statistics                    #----
#---------------------------------------------------------------------#

#10. Calculate r-squared and F-statistics

fstat_r2_func <- function(exposure_name) {
  total_r2 <- sum(dat[dat$exposure==exposure_name,]$rsq.exposure)
  totalSNPs <- nrow(dat[dat$exposure==exposure_name,])
   min_F <- min(dat[dat$exposure==exposure_name,]$F, na.rm = T)
   max_F <- max(dat[dat$exposure==exposure_name,]$F, na.rm = T)
   mean_F <- mean(dat[dat$exposure==exposure_name,]$F, na.rm = T)
   effectiveSNPs <- nrow(dat[dat$exposure==exposure_name & dat$mr_keep==T,])
   total_r2_effective <- sum(dat[dat$exposure==exposure_name & dat$mr_keep==T,]$rsq.exposure)
   min_F_effective <- min(dat[dat$exposure==exposure_name & dat$mr_keep==T,]$F, na.rm = T)
   max_F_effective <- max(dat[dat$exposure==exposure_name & dat$mr_keep==T,]$F, na.rm = T)
   mean_F_effective <- mean(dat[dat$exposure==exposure_name & dat$mr_keep==T,]$F, na.rm = T)
   
   paste0("for ", exposure_name, " with N SNPS:", totalSNPs, ", the total r2=", total_r2, ", the mean F-statistic=", 
          mean_F, ", the min F=", min_F, ", and the max F=", max_F, ". EFFECTIVE SNPS:", effectiveSNPs, ", the total effective r2=", 
          total_r2_effective, ", the mean effective F-statistic=", mean_F_effective, ", the min effective F=", 
          min_F_effective, ", and the max effective F=", max_F_effective)
   
}

BMI_strict <- fstat_r2_func("BMI_strict")
WHR_strict <- fstat_r2_func("WHR_strict")
WC <- fstat_r2_func("WC")

adipose <- fstat_r2_func("adipose-BMI")
brain <- fstat_r2_func("brain-BMI")
adult <- fstat_r2_func("adult BMI")
child <- fstat_r2_func("child BMI")
BF <- fstat_r2_func("Body fat percentage")
FA <- fstat_r2_func("FA")
UFA <- fstat_r2_func("UFA")
PC1 <- fstat_r2_func("PC1")
PC2 <- fstat_r2_func("PC2")
PC3 <- fstat_r2_func("PC3")
PC4 <- fstat_r2_func("PC4")

all <- c(BMI_strict, WHR_strict, WC, adipose, brain, adult, child, BF,
         FA, UFA, PC1, PC2, PC3, PC4)
all

#"for BMI_strict with N SNPS:451, the total r2=0.050157172518088, the mean F-statistic=79.6096490037, the min F=32.801652892562, and the max F=1270.71280276817. EFFECTIVE SNPS:442, the total effective r2=0.0479821345449918, the mean effective F-statistic=77.3903696745938, the min effective F=32.801652892562, and the max effective F=844.05540166205"       
#"for WHR_strict with N SNPS:276, the total r2=0.0319247572089752, the mean F-statistic=72.7562931036939, the min F=32.9113573407202, and the max F=820.041322314049. EFFECTIVE SNPS:267, the total effective r2=0.0310638258337821, the mean effective F-statistic=73.2356937451418, the min effective F=32.9113573407202, and the max effective F=820.041322314049"
#"for WC with N SNPS:38, the total r2=0.00990844253511101, the mean F-statistic=60.5557841679886, the min F=29.3402777777778, and the max F=447.020408163265. EFFECTIVE SNPS:37, the total effective r2=0.0079861546052696, the mean effective F-statistic=50.1107943302784, the min effective F=29.3402777777778, and the max effective F=144"
###for WC with N SNPS:368, the total r2=0.0455859755702108, the mean F-statistic=57.2656885229747, the min F=29.7606843509764, and the max F=940.083866154503. EFFECTIVE SNPS:353, the total effective r2=0.0442269385855828, the mean effective F-statistic=57.919583716967, the min effective F=29.7606843509764, and the max effective F=940.083866154503
#"for adipose-BMI with N SNPS:84, the total r2=0.00763666488502167, the mean F-statistic=61.9452687575298, the min F=29.641975308642, and the max F=270.19140625. EFFECTIVE SNPS:81, the total effective r2=0.00745644124381982, the mean effective F-statistic=62.7236252876156, the min effective F=29.641975308642, and the max effective F=270.19140625"                  
#"for brain-BMI with N SNPS:136, the total r2=0.0120502032718485, the mean F-statistic=60.3718144227282, the min F=29.369406867846, and the max F=270.19140625. EFFECTIVE SNPS:133, the total effective r2=0.0118635198393921, the mean effective F-statistic=60.7772665418973, the min effective F=29.369406867846, and the max effective F=270.19140625"                    
#"for adult BMI with N SNPS:335, the total r2=0.0435934497838226, the mean F-statistic=58.989162193629, the min F=29.8818840927864, and the max F=1108.84561855914. EFFECTIVE SNPS:324, the total effective r2=0.0422049441685177, the mean effective F-statistic=59.0494481921754, the min effective F=29.8818840927864, and the max effective F=1108.84561855914"
#"for child BMI with N SNPS:203, the total r2=0.034325733055672, the mean F-statistic=76.6723652668551, the min F=27.8902226522023, and the max F=1102.20315561404. EFFECTIVE SNPS:198, the total effective r2=0.0338445501448658, the mean effective F-statistic=77.5071251681359, the min effective F=28.4262399131224, and the max effective F=1102.20315561404"           
#"for Body fat percentage with N SNPS:371, the total r2=0.0478272436001003, the mean F-statistic=58.6209622508114, the min F=29.8030765762561, and the max F=681.934702635416. EFFECTIVE SNPS:360, the total effective r2=0.0465727885781771, the mean effective F-statistic=58.8277357402828, the min effective F=29.8030765762561, and the max effective F=681.934702635416"
#"for FA with N SNPS:32, the total r2=0.00462913251951022, the mean F-statistic=65.2777777777778, the min F=25, and the max F=400. EFFECTIVE SNPS:31, the total effective r2=0.00440749984021492, the mean effective F-statistic=64.15770609319, the min effective F=25, and the max effective F=400"                                                                         
#"for UFA with N SNPS:28, the total r2=0.00808758622923102, the mean F-statistic=130.357142857143, the min F=25, and the max F=400. EFFECTIVE SNPS:27, the total effective r2=0.00786595354993573, the mean effective F-statistic=131.481481481481, the min effective F=25, and the max effective F=400"                                                                      
#"for PC1 with N SNPS:28, the total r2=0.159943295080727, the mean F-statistic=54.4890371705766, the min F=27.5625, and the max F=302.457466918715. EFFECTIVE SNPS:28, the total effective r2=0.159943295080727, the mean effective F-statistic=54.4890371705766, the min effective F=27.5625, and the max effective F=302.457466918715"                                      
#"for PC2 with N SNPS:84, the total r2=0.03554990878132, the mean F-statistic=53.7693844295472, the min F=29.7520661157025, and the max F=210.534409842368. EFFECTIVE SNPS:81, the total effective r2=0.0343758484855078, the mean effective F-statistic=53.8806087673288, the min effective F=29.7520661157025, and the max effective F=210.534409842368"                    
#"for PC3 with N SNPS:27, the total r2=0.00853437109431019, the mean F-statistic=41.2126323909828, the min F=30.25, and the max F=82.2606814494321. EFFECTIVE SNPS:27, the total effective r2=0.00853437109431019, the mean effective F-statistic=41.2126323909828, the min effective F=30.25, and the max effective F=82.2606814494321"                                      
#"for PC4 with N SNPS:10, the total r2=0.246782104564571, the mean F-statistic=42.1517690782406, the min F=30.25, and the max F=97.6608996539792. EFFECTIVE SNPS:10, the total effective r2=0.246782104564571, the mean effective F-statistic=42.1517690782406, the min effective F=30.25, and the max effective F=97.6608996539792"

###----------------------------------------------

  
  
  
  
  
  
  
  
  
  
  
  
#---------------------------------------------------------------------#
#                R2 and F-statistic (other methods) using clumped data instead of harmonised data                   #----
#---------------------------------------------------------------------#
#load data
dat_BMI_strict <- fread("results/clumped_exposure_files/clumped_BMI.csv", header = T) %>% filter(., pval.exposure<5e-9)
dat_WHR_strict <- fread("results/clumped_exposure_files/clumped_WHR.csv", header = T) %>% filter(., pval.exposure<5e-9)
dat_WC <- fread("results/clumped_exposure_files/clumped_WC.csv", header = T) %>% mutate(samplesize.exposure=232101)


dat_PC1_strict <- fread("results/clumped_exposure_files/clumped_PC1.csv", header = T) 
dat_PC2_strict <- fread("results/clumped_exposure_files/clumped_PC2.csv", header = T) 
dat_PC3_strict <- fread("results/clumped_exposure_files/clumped_PC3.csv", header = T) 
dat_PC4_strict <- fread("results/clumped_exposure_files/clumped_PC4.csv", header = T) 
dat_FA_strict <- fread("results/clumped_exposure_files/clumped_FA.csv", header = T) %>% mutate(samplesize.exposure=451099)
dat_UFA_strict <- fread("results/clumped_exposure_files/clumped_UFA.csv", header = T) %>% mutate(samplesize.exposure=451099)

dat_richardson <- extract_instruments(
  c("ieu-b-5107", "ukb-b-2303"),
  p1 = 5e-08,
  clump = TRUE,
  p2 = 5e-08,
  r2 = 0.001,
  kb = 10000,
  access_token = ieugwasr::check_access_token(),
  force_server = FALSE
)

dat_richardson <- clump_data(dat_richardson, clump_p1 = 5e-8, clump_p2 = 5e-8)
dat_child_richardson <- dat_richardson %>% filter(., exposure=="Comparative body size at age 10, Males and Females || id:ieu-b-5107") %>% mutate(samplesize.exposure=453169)
dat_adult_richardson <-dat_richardson %>% filter(., exposure=="Body mass index (BMI) || id:ukb-b-2303") 


dat_fat <- extract_instruments(
  c("ukb-b-8909"),
  p1 = 5e-08,
  clump = TRUE,
  p2 = 5e-08,
  r2 = 0.001,
  kb = 10000,
  access_token = ieugwasr::check_access_token(),
  force_server = FALSE
)

dat_fat <- clump_data(dat_fat, clump_p1 = 5e-8, clump_p2 = 5e-8)


brain_dat <- read.xlsx("data/GIANT_Pulit2018/tissue_partitioned_BMI.xlsx") %>% filter(., Tissue=="brain" | Tissue=="both") %>% format_data(type = "exposure") %>% mutate(samplesize.exposure=681275)
adipose_dat <- read.xlsx("data/GIANT_Pulit2018/tissue_partitioned_BMI.xlsx") %>% filter(., Tissue=="adipose" | Tissue=="both") %>% format_data(type = "exposure") %>% mutate(samplesize.exposure=681275)


#create function

fstat_r2_func2 <- function(exposure_name, exp_dat) {
  exp_dat$r2 <- (2 * (exp_dat$beta.exposure^2) * exp_dat$eaf.exposure * (1 - exp_dat$eaf.exposure)) /
    (2 * (exp_dat$beta.exposure^2) * exp_dat$eaf.exposure * (1 - exp_dat$eaf.exposure) +
       2 * exp_dat$samplesize.exposure * exp_dat$eaf.exposure * 
       (1 - exp_dat$eaf.exposure) * exp_dat$se.exposure^2)
  exp_dat$F <- exp_dat$r2 * (exp_dat$samplesize.exposure - 2) / (1 - exp_dat$r2)
  
  total_r2 <- sum(exp_dat$r2, na.rm = T)
  min_r2 <- min(exp_dat$r2, na.rm = T)
  max_r2 <- max(exp_dat$r2, na.rm = T)
  min_F <- min(exp_dat$F, na.rm = T)
  max_F <- max(exp_dat$F, na.rm = T)
  mean_F <- mean(exp_dat$F, na.rm = T)

  paste0("for ", exposure_name, ", the total r2=", total_r2," (range ", min_r2, "-",max_r2, ") and the mean F-statistic=", mean_F, " (range", min_F, "-",max_F, ")")
  
}

BMI_strict <- fstat_r2_func2("BMI_strict", dat_BMI_strict)
WHR_strict <- fstat_r2_func2("WHR_strict", dat_WHR_strict)
WC <- fstat_r2_func2("WC", dat_WC)

height_strict <- fstat_r2_func2("height_strict", dat_height_strict)
PC1_strict <- fstat_r2_func2("PC1_strict", dat_PC1_strict)
PC2_strict <- fstat_r2_func2("PC2_strict", dat_PC2_strict)
PC3_strict <- fstat_r2_func2("PC3_strict", dat_PC3_strict)
PC4_strict <- fstat_r2_func2("PC4_strict", dat_PC4_strict)
FA_strict <- fstat_r2_func2("FA_strict", dat_FA_strict)
UFA_strict <- fstat_r2_func2("UFA_strict", dat_UFA_strict)

child_strict <- fstat_r2_func2("child_strict", dat_child_richardson)
adult_strict <- fstat_r2_func2("adult_strict", dat_adult_richardson)

fat_strict <- fstat_r2_func2("fat_strict", dat_fat)


brain_strict <- fstat_r2_func2("brain_strict", brain_dat)
adipose_strict <- fstat_r2_func2("adipose_strict", adipose_dat)


#we need to multiple r2 by 100 to obtain %
