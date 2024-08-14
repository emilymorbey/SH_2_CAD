library(tidyverse)
library(readxl)
library(MendelianRandomization)

setwd("C:/Users/emorb/OneDrive - University of Cambridge/PhD/MR/Testosterone_CAD_MR/Testosterone CAD MR R files/TestosteroneCAD")

##################################################################################

# HARMONISATION AND MR

##################################################################################

# looking at the allele matching and frequencies etc.
M_T_proxies_output <- read_excel("not found inputs/SNPs_M_Testosterone_AND_CRP.xlsx", sheet = "T&P E&O")
View(M_T_proxies_output)




allele_matching <- select(M_T_proxies_output, "CHR...2", "BP...6", "SNP", "ALLELE1", "ALLELE0", "BETA", "SE", "effect_allele", "other_allele", "beta", "standard_error" )

# renaming the columns for ease of use 
allele_matching <- allele_matching %>%
  rename(
    SNP_T = "SNP",
    other_allele_T = "ALLELE1",
    reference_allele_T = "ALLELE0",
    BETA_T = "BETA",
    SE_T = "SE",
    reference_allele_CRP = "effect_allele",
    other_allele_CRP = "other_allele",
    beta_CRP = "beta",
    se_CRP = "standard_error"
  )

# Removing SNP that was pleitropic in CAD MR 
allele_matching <- allele_matching[!allele_matching$SNP_T == "rs56196860", ]
M_T_proxies_output <- M_T_proxies_output[!M_T_proxies_output$TARGET == "rs56196860", ]


# identify trait increasing allele for SHBG

allele_matching$T_inc_allele <- if_else(allele_matching$BETA_T<0, allele_matching$reference_allele_T, 
                                        allele_matching$other_allele_T)

allele_matching$BETA_T <- as.numeric(allele_matching$BETA_T)
allele_matching$ABS_BETA_T <- abs(allele_matching$BETA_T)

# harmonising so the effect alleles for CAD and CRP are the same
# changing the betas here 

allele_matching$beta_CRP <- as.numeric(allele_matching$beta_CRP)

# CAPITALISE ALL THE ALLELE VALUES FOR DBP

allele_matching$other_allele_CRP <- toupper(allele_matching$other_allele_CRP)
allele_matching$reference_allele_CRP <- toupper(allele_matching$reference_allele_CRP)

allele_matching$HARM_BETA_CRP <- if_else(allele_matching$T_inc_allele!=allele_matching$other_allele_CRP,
                                         allele_matching$beta_CRP*-1, allele_matching$beta_CRP)




plot(allele_matching$ABS_BETA_T, allele_matching$HARM_BETA_CRP)


M_T_proxies_output$StdError <- as.numeric(M_T_proxies_output$standard_error)
IVW_weights <- M_T_proxies_output$standard_error^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_BETA_CRP ~ allele_matching$ABS_BETA_T- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red")
summary_model <- summary(inverse_weighted_LR)
summary_model



allele_matching$ABS_BETA_T <- as.numeric(allele_matching$ABS_BETA_T)
allele_matching$se_CRP <- as.numeric(allele_matching$se_CRP)
allele_matching$SE_T <- as.numeric(allele_matching$SE_T)



MRObject = mr_input(bx = allele_matching$ABS_BETA_T, bxse = allele_matching$SE_T, 
                    by = allele_matching$HARM_BETA_CRP, byse = allele_matching$se_CRP, snps = allele_matching$SNP_T)

mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)

mr_allmethods(MRObject)
mr_allmethods(MRObject, method="main")


plot(allele_matching$ABS_BETA_T, allele_matching$HARM_BETA_CRP  ,
     xlab = "SNP effect on testosterone",  # Replace with your desired x-axis label
     ylab = "SNP effect on CRP",
     main = "Male Testosterone")  # Replace with your desired y-axis label

M_T_proxies_output$standard_error <- as.numeric(M_T_proxies_output$standard_error)
IVW_weights <- M_T_proxies_output$standard_error^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_BETA_CRP ~ allele_matching$ABS_BETA_T - 1, weights = IVW_weights)
summary(inverse_weighted_LR)

abline(inverse_weighted_LR, col = "red")



##################################################################################




#######################################################################################

### other interprative graphs #############

mr_plot(MRObject, interactive=FALSE, labels=TRUE)
mr_forest(MRObject, ordered=TRUE)
mr_loo(MRObject)
mr_funnel(MRObject)
