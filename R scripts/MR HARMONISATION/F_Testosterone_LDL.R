library(tidyverse)
library(readxl)
library(MendelianRandomization)

##################################################################################

# HARMONISATION AND MR

##################################################################################

# looking at the allele matching and frequencies etc.
F_T_proxies_output <- read_excel("not found inputs/SNPs_F_Testosterone_AND_LDL.xlsx", sheet = "T&P E&O")



allele_matching <- select(F_T_proxies_output, "SNP.x", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE.x", "REF", "ALT", "POOLED_ALT_AF", "EFFECT_SIZE", "SE.y" )

# renaming the columns for ease of use 
allele_matching <- allele_matching %>%
  rename(
    SNP_T = "SNP.x",
    other_allele_T = "ALLELE1",
    reference_allele_T = "ALLELE0",
    A1FREQ_T = "A1FREQ",
    BETA_T = "BETA",
    SE_T = "SE.x",
    reference_allele_LDL = "REF",
    other_allele_LDL = "ALT",
    eaf_LDL = "POOLED_ALT_AF",
    beta_LDL = "EFFECT_SIZE",
    se_LDL = "SE.y"
  )

# identify trait increasing allele for SHBG

allele_matching$T_inc_allele <- if_else(allele_matching$BETA_T<0, allele_matching$reference_allele_T, 
                                        allele_matching$other_allele_T)

allele_matching$BETA_T <- as.numeric(allele_matching$BETA_T)
allele_matching$ABS_BETA_T <- abs(allele_matching$BETA_T)

# harmonising so the effect alleles for CAD and SHBG are the same
# changing the betas here 

allele_matching$beta_LDL <- as.numeric(allele_matching$beta_LDL)

allele_matching$HARM_BETA_LDL <- if_else(allele_matching$T_inc_allele!=allele_matching$other_allele_LDL,
                                         allele_matching$beta_LDL*-1, allele_matching$beta_LDL)




plot(allele_matching$ABS_BETA_T, allele_matching$HARM_BETA_LDL)


F_T_proxies_output$SE.y <- as.numeric(F_T_proxies_output$SE.y)
IVW_weights <- F_T_proxies_output$SE.y^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_BETA_LDL ~ allele_matching$ABS_BETA_T- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red")
summary_model <- summary(inverse_weighted_LR)
summary_model



allele_matching$ABS_BETA_T <- as.numeric(allele_matching$ABS_BETA_T)
allele_matching$se_LDL <- as.numeric(allele_matching$se_LDL)
allele_matching$SE_T <- as.numeric(allele_matching$SE_T)



MRObject = mr_input(bx = allele_matching$ABS_BETA_T, bxse = allele_matching$SE_T, 
                    by = allele_matching$HARM_BETA_LDL, byse = allele_matching$se_LDL, snps = allele_matching$SNP_T)

mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)



mr_plot(mr_allmethods(MRObject))



plot(allele_matching$ABS_BETA_T, allele_matching$HARM_BETA_LDL  ,
     xlab = "SNP effect on Testosterone",  # Replace with your desired x-axis label
     ylab = "SNP effect on LDL-c",
     main = "Female Testosterone")  # Replace with your desired y-axis label

F_T_proxies_output$SE.y <- as.numeric(F_T_proxies_output$SE.y)
IVW_weights <- F_T_proxies_output$SE.y^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_BETA_LDL ~ allele_matching$ABS_BETA_T - 1, weights = IVW_weights)
summary(inverse_weighted_LR)

abline(inverse_weighted_LR, col = "red")


#######################################################################################

### other interprative graphs #############

mr_plot(MRObject, interactive=FALSE, labels=TRUE)
mr_forest(MRObject, ordered=TRUE)
mr_loo(MRObject)
mr_funnel(MRObject)
















