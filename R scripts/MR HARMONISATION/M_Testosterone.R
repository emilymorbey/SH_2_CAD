
install.packages("xfun")
##################################################################################

# HARMONISATION AND MR

##################################################################################
library(tidyverse)
library(readxl)
library(MendelianRandomization)

# looking at the allele matching and frequencies etc.
M_T_proxies_output <- read_excel("not found inputs/SNPs_M_Testosterone_AND_CAD.xlsx", sheet = "FREE T E&O")



allele_matching <- select(M_T_proxies_output, "SNP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "reference_allele", "other_allele", "eaf", "male_beta", "male_se" )

# renaming the columns for ease of use 
allele_matching <- allele_matching %>%
  rename(
    SNP_T = "SNP",
    Effect_allele_T = "ALLELE1",
    Reference_allele_T = "ALLELE0",
    EA_FREQ_T = "A1FREQ",
    BETA_T = "BETA",
    SE_T = "SE",
    Reference_allele_CAD = "reference_allele",
    Effect_allele_CAD = "other_allele",
    EA_FREQ_CAD = "eaf",
    male_beta_CAD = "male_beta",
    male_se_CAD = "male_se"
  )

# identify trait increasing allele for SHBG

allele_matching$T_inc_allele <- if_else(allele_matching$BETA_T<0, allele_matching$Reference_allele_T, 
                                           allele_matching$Effect_allele_T)

allele_matching$BETA_T <- as.numeric(allele_matching$BETA_T)
allele_matching$ABS_BETA_T <- abs(allele_matching$BETA_T)

# harmonising so the effect alleles for CAD and SHBG are the same
# changing the betas here 

allele_matching$male_beta_CAD <- as.numeric(allele_matching$male_beta_CAD)

allele_matching$HARM_MALE_BETA_CAD <- if_else(allele_matching$T_inc_allele!=allele_matching$Effect_allele_CAD,
                                              allele_matching$male_beta_CAD*-1, allele_matching$male_beta_CAD)




plot(allele_matching$ABS_BETA_T, allele_matching$HARM_MALE_BETA_CAD)
allele_matching$male_se_CAD <- as.numeric(allele_matching$male_se_CAD)
IVW_weights <- allele_matching$male_se_CAD^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_T- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red")
summary_model <- summary(inverse_weighted_LR)
summary_model



M_T_proxies_output$male_beta <- as.numeric(M_T_proxies_output$male_beta)

allele_matching$ABS_BETA_T <- as.numeric(allele_matching$ABS_BETA_T)
allele_matching$male_se_CAD <- as.numeric(allele_matching$male_se_CAD)
allele_matching$SE_T <- as.numeric(allele_matching$SE_T)

MRObject = mr_input(bx = allele_matching$ABS_BETA_T, bxse = allele_matching$SE_T, 
                    by = allele_matching$HARM_MALE_BETA_CAD, byse = allele_matching$male_se_CAD, snps = allele_matching$SNP_T)

mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)

mr_allmethods(MRObject)
mr_plot(mr_allmethods(MRObject))



plot(allele_matching$ABS_BETA_T, allele_matching$HARM_MALE_BETA_CAD,
     xlab = "SNP effect on Testosterone",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Male Testosterone")  # Replace with your desired y-axis label

M_T_proxies_output$male_se <- as.numeric(M_T_proxies_output$male_se)
IVW_weights <- M_T_proxies_output$male_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_T - 1, weights = IVW_weights)
summary(inverse_weighted_LR)

abline(inverse_weighted_LR, col = "red")





#######################################################################################

### other interprative graphs #############

mr_plot(MRObject, interactive=TRUE, labels=TRUE)
mr_forest(MRObject, ordered=TRUE)
mr_loo(MRObject)
mr_funnel(MRObject)
??mr_funnel


mr_funnel(MRObject) + ggtitle("Funnel plot for male testosterone and CAD")


################ RUNNING WITHOUT OUTLIER ############################



# leaving out the SNP that is causing the problem

allele_matching <- allele_matching[!allele_matching$SNP_T == "rs56196860", ]
M_T_proxies_output <- M_T_proxies_output[!M_T_proxies_output$Target == "rs56196860", ]
# identify trait increasing allele for SHBG

allele_matching$T_inc_allele <- if_else(allele_matching$BETA_T<0, allele_matching$Reference_allele_T, 
                                        allele_matching$Effect_allele_T)

allele_matching$BETA_T <- as.numeric(allele_matching$BETA_T)
allele_matching$ABS_BETA_T <- abs(allele_matching$BETA_T)

# harmonising so the effect alleles for CAD and SHBG are the same
# changing the betas here 

allele_matching$male_beta_CAD <- as.numeric(allele_matching$male_beta_CAD)

allele_matching$HARM_MALE_BETA_CAD <- if_else(allele_matching$T_inc_allele!=allele_matching$Effect_allele_CAD,
                                              allele_matching$male_beta_CAD*-1, allele_matching$male_beta_CAD)




plot(allele_matching$ABS_BETA_T, allele_matching$HARM_MALE_BETA_CAD)
M_T_proxies_output$male_se <- as.numeric(M_T_proxies_output$male_se)
IVW_weights <- M_T_proxies_output$male_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_T- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red")
summary_model <- summary(inverse_weighted_LR)
summary_model



M_T_proxies_output$male_beta <- as.numeric(M_T_proxies_output$male_beta)

allele_matching$ABS_BETA_T <- as.numeric(allele_matching$ABS_BETA_T)
allele_matching$male_se_CAD <- as.numeric(allele_matching$male_se_CAD)
allele_matching$SE_T <- as.numeric(allele_matching$SE_T)

MRObject = mr_input(bx = allele_matching$ABS_BETA_T, bxse = allele_matching$SE_T, 
                    by = allele_matching$HARM_MALE_BETA_CAD, byse = allele_matching$male_se_CAD, snps = allele_matching$SNP_T)

mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)

mr_allmethods(MRObject)
mr_plot(mr_allmethods(MRObject))
mr_plot(MRObject, interactive=TRUE, labels=TRUE)


plot(allele_matching$ABS_BETA_T, allele_matching$HARM_MALE_BETA_CAD,
     xlab = "SNP effect on Testosterone",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Male Testosterone")  # Replace with your desired y-axis label

M_T_proxies_output$male_se <- as.numeric(M_T_proxies_output$male_se)
IVW_weights <- M_T_proxies_output$male_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_T - 1, weights = IVW_weights)
summary(inverse_weighted_LR)

abline(inverse_weighted_LR, col = "red")
