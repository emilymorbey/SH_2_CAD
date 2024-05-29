


################################################################################

# HARMONISATION AND MR

##################################################################################
library(tidyverse)
library(readxl)
library(MendelianRandomization)

# looking at the allele matching and frequencies etc.
F_SHBG_proxies_output <- read_excel("not found inputs/SNPs_F_SHBG_AND_CAD.xlsx", sheet = "T&P E&O")
F_SHBG_proxies_output <- F_SHBG_proxies_output[-1,]


allele_matching <- select(F_SHBG_proxies_output, "SNP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "reference_allele", "other_allele", "eaf", "female_beta", "female_se" )

# renaming the columns for ease of use 
allele_matching <- allele_matching %>%
  rename(
    SNP_SHBG = "SNP",
    ALLELE1_SHBG = "ALLELE1",
    ALLELE0_SHBG = "ALLELE0",
    A1FREQ_SHBG = "A1FREQ",
    BETA_SHBG = "BETA",
    SE_SHBG = "SE",
    reference_allele_CAD = "reference_allele",
    other_allele_CAD = "other_allele",
    eaf_CAD = "eaf",
    female_beta_CAD = "female_beta",
    female_se_CAD = "female_se"
  )

# identify trait increasing allele for SHBG

allele_matching$SHBG_inc_allele <- if_else(allele_matching$BETA_SHBG<0, allele_matching$ALLELE0_SHBG, 
                                           allele_matching$ALLELE1_SHBG)

allele_matching$BETA_SHBG <- as.numeric(allele_matching$BETA_SHBG)
allele_matching$ABS_BETA_SHBG <- abs(allele_matching$BETA_SHBG)

# harmonising so the effect alleles for CAD and SHBG are the same
# changing the betas here 

allele_matching$female_beta_CAD <- as.numeric(allele_matching$female_beta_CAD)

allele_matching$HARM_FEMALE_BETA_CAD <- if_else(allele_matching$SHBG_inc_allele!=allele_matching$other_allele_CAD,
                                              allele_matching$female_beta_CAD*-1, allele_matching$female_beta_CAD)




plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_FEMALE_BETA_CAD)
F_SHBG_proxies_output$female_se <- as.numeric(F_SHBG_proxies_output$female_se)
IVW_weights <- F_SHBG_proxies_output$female_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_FEMALE_BETA_CAD ~ allele_matching$ABS_BETA_SHBG- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red")
summary_model <- summary(inverse_weighted_LR)
summary_model



F_SHBG_proxies_output$female_beta <- as.numeric(F_SHBG_proxies_output$female_beta)

allele_matching$ABS_BETA_SHBG <- as.numeric(allele_matching$ABS_BETA_SHBG)
allele_matching$female_se_CAD <- as.numeric(allele_matching$female_se_CAD)
allele_matching$SE_SHBG <- as.numeric(allele_matching$SE_SHBG)

MRObject = mr_input(bx = allele_matching$ABS_BETA_SHBG, bxse = allele_matching$SE_SHBG, 
                    by = allele_matching$HARM_FEMALE_BETA_CAD, byse = allele_matching$female_se_CAD, snps = allele_matching$SNP_SHBG)

mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)


mr_allmethods(MRObject)
mr_plot(MRObject, error = TRUE, line = "allmethods", interactive = FALSE, labels=FALSE)

mr_plot(mr_allmethods(MRObject))





plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_FEMALE_BETA_CAD,
     xlab = "SNP effect on SHBG",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Female SHBG")  # Replace with your desired y-axis label

F_SHBG_proxies_output$female_se <- as.numeric(F_SHBG_proxies_output$female_se)
IVW_weights <- F_SHBG_proxies_output$female_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_FEMALE_BETA_CAD ~ allele_matching$ABS_BETA_SHBG - 1, weights = IVW_weights)
summary(inverse_weighted_LR)

abline(inverse_weighted_LR, col = "red")





# forest plot

mr_forest(MRObject, alpha = 0.05, snp_estimates = FALSE, methods = "ivw", ordered = FALSE)
mr_plot(MRObject, interactive=FALSE, labels=TRUE)
  
allele_matching %>%
  filter(SNP_SHBG=="rs6258")
