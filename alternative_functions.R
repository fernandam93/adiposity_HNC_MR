#2. Extract SNP list for BMI from HEADSpAcE GWAS summary data 

#if SNPs have been extracted beforehand using linux, this is the best alternative
outcome_harmonisation_func_alternative <- function(outcome_file, outcome_name, exposure_dat, exposure_name,
                                                   snp="V19", beta="BETA", se="SE",
                                                   eaf="A1_FREQ", effect_allele="A1",other_allele="AX",
                                                   pval="P", chr="CHR", pos="BP", sep = "\t") {
  
  outcome_stats <- read_outcome_data(snps = exposure_dat$SNP,
                                     filename = outcome_file,
                                     sep = sep,
                                     snp_col = snp,
                                     beta_col = beta,
                                     se_col = se,
                                     effect_allele_col = effect_allele,
                                     other_allele_col = other_allele,
                                     eaf_col = eaf,
                                     pval_col = pval,
                                     min_pval = 1e-200,
                                     log_pval = FALSE,
                                     chr_col = chr,
                                     pos_col = pos,
                                     samplesize_col = "N")
  outcome_stats$outcome <- outcome_name
  
  #harmonisation
  dat <- harmonise_data(exposure_dat, outcome_stats)
  dat$exposure <- exposure_name
  
  #save
  write.csv(dat, paste0('results/harmonised_data_files/harmonised_', exposure_name, '_', outcome_name,'.csv'), row.names = F)
  
  dat
  
}


#if the outcome is on MR-BASE, this is the best alternative
outcome_harmonisation_func_mrbase <- function(outcome_mrbase, outcome_name, exposure_dat, exposure_name) {
  
  outcome_stats <- extract_outcome_data(
                                    snps = exposure_dat$SNP,
                                    outcomes = outcome_mrbase,
                                    opengwas_jwt = ieugwasr::get_opengwas_jwt()
                                  )
  outcome_stats$outcome <- outcome_name
  
  #harmonisation
  dat <- harmonise_data(exposure_dat, outcome_stats)
  dat$exposure <- exposure_name
  
  #save
  write.csv(dat, paste0('results/harmonised_data_files/harmonised_', exposure_name, '_', outcome_name,'.csv'), row.names = F)
  
  dat
  
}

