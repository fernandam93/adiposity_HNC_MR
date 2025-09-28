# adiposity_HNC_MR
This code relates to the project titled "Reassessing the link between adiposity and head and neck cancer: a Mendelian randomization study", published in eLife (DOI: https://doi.org/10.7554/eLife.106075).

## Step 1: Create exposure datasets

1.headspace_MR_exposure_prep.R

## Step 2: Harmonise exposure-outcome datasets and run two-sample MR analysis

2.harmonisation_MR_noukb_inc_proxies.R

## Step 3: Create main plots

3.headspace_MR_plots.R

## Step 4: Secondary analyses (heterogeneity across HNC subsites for main analyses, and MR using secondary adiposity traits as exposures)

4.0.bmi_subtype_heterogeneity.R

4.1.bodyshape_MR.R

4.2.bodysize_MR.R

4.3.partitionedBMI_MR.R

4.4.UFA_FA_MR.R

## Step 5: MVMR analysis accounting for smoking behaviour

5.1.bmi_smoking_mvmr.R

5.2.bmi_csi_mvmr_plot.R

5.3.bmi_si_mvmr_plot.R

5.4.si_csi_univariable_IVW.R

## Step 6: MR-PRESSO

6.MRPRESSO_noukb_inc_proxies.R

## Step 7: MR using Steiger filtering to remove SNPs that are closer to heritable confounders than BMI

7.steiger_filtering_MR.R

## Step 8: CAUSE 

8.CAUSE_BMI_HNC.R

## Step 9: MR-Clust

9.bmi_HNC_mrclust.R

## Additional scripts
**R2 and F-statistics were calculated for the instruments used in the analyses**

adiposity_r2_fstats.R

**File containing functions used across scripts**

alternative_functions.R
