library(tidyverse)
library(readxl)
library(MendelianRandomization)


##################################################################################

# HARMONISATION AND MR

##################################################################################

# looking at the allele matching and frequencies etc.
F_T_proxies_output <- read_excel("not found inputs/SNPs_F_Testosterone_AND_CAD.xlsx", sheet = "T&P E&O")



allele_matching <- select(F_T_proxies_output, "SNP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "reference_allele", "other_allele", "eaf", "female_beta", "female_se" )

# renaming the columns for ease of use 
allele_matching <- allele_matching %>%
  rename(
    SNP_T = "SNP",
    ALLELE1_T = "ALLELE1",
    ALLELE0_T = "ALLELE0",
    A1FREQ_T = "A1FREQ",
    BETA_T = "BETA",
    SE_T = "SE",
    reference_allele_CAD = "reference_allele",
    other_allele_CAD = "other_allele",
    eaf_CAD = "eaf",
    female_beta_CAD = "female_beta",
    female_se_CAD = "female_se"
  )

# identify trait increasing allele for SHBG

allele_matching$T_inc_allele <- if_else(allele_matching$BETA_T<0, allele_matching$ALLELE0_T, 
                                        allele_matching$ALLELE1_T)

allele_matching$BETA_T <- as.numeric(allele_matching$BETA_T)
allele_matching$ABS_BETA_T <- abs(allele_matching$BETA_T)

# harmonising so the effect alleles for CAD and SHBG are the same
# changing the betas here 

allele_matching$female_beta_CAD <- as.numeric(allele_matching$female_beta_CAD)

allele_matching$HARM_FEMALE_BETA_CAD <- if_else(allele_matching$T_inc_allele!=allele_matching$other_allele_CAD,
                                              allele_matching$female_beta_CAD*-1, allele_matching$female_beta_CAD)




plot(allele_matching$ABS_BETA_T, allele_matching$HARM_FEMALE_BETA_CAD)
F_T_proxies_output$female_se <- as.numeric(F_T_proxies_output$female_se)
IVW_weights <- F_T_proxies_output$female_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_FEMALE_BETA_CAD ~ allele_matching$ABS_BETA_T- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red")
summary_model <- summary(inverse_weighted_LR)
summary_model



F_T_proxies_output$female_beta <- as.numeric(F_T_proxies_output$female_beta)

allele_matching$ABS_BETA_T <- as.numeric(allele_matching$ABS_BETA_T)
allele_matching$female_se_CAD <- as.numeric(allele_matching$female_se_CAD)
allele_matching$SE_T <- as.numeric(allele_matching$SE_T)



MRObject = mr_input(bx = allele_matching$ABS_BETA_T, bxse = allele_matching$SE_T, 
                    by = allele_matching$HARM_FEMALE_BETA_CAD, byse = allele_matching$female_se_CAD, snps = allele_matching$SNP_T)

mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)


mr_allmethods(MRObject)
mr_plot(mr_allmethods(MRObject))



plot(allele_matching$ABS_BETA_T, allele_matching$HARM_FEMALE_BETA_CAD,
     xlab = "SNP effect on Testosterone",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Female Testosterone")  # Replace with your desired y-axis label

F_T_proxies_output$female_se <- as.numeric(F_T_proxies_output$female_se)
IVW_weights <- F_T_proxies_output$female_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_FEMALE_BETA_CAD ~ allele_matching$ABS_BETA_T - 1, weights = IVW_weights)
summary(inverse_weighted_LR)

abline(inverse_weighted_LR, col = "red")





##################################################################################




#######################################################################################

### other interprative graphs #############

mr_plot(MRObject, interactive=FALSE, labels=TRUE)
mr_forest(MRObject, ordered=TRUE)
mr_loo(MRObject)
mr_funnel(MRObject)









duplicates <- allele_matching$SNP_T[duplicated(allele_matching$SNP_T)]
print(duplicates)





