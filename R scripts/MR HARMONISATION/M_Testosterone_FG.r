##################################################################################

# HARMONISATION AND MR

##################################################################################
library(tidyverse)
library(readxl)
library(MendelianRandomization)


setwd("C:/Users/emorb/OneDrive - University of Cambridge/PhD/MR/Testosterone_CAD_MR/Testosterone CAD MR R files")

# looking at the allele matching and frequencies etc.
M_T_proxies_output <- read_excel("TestosteroneCAD/not found inputs/SNPs_M_T_FG.xlsx", sheet = "T&P E&O2")
View(M_T_proxies_output)

M_T_proxies_output <- M_T_proxies_output[!M_T_proxies_output$rsid == "rs56196860", ]

allele_matching <- select(M_T_proxies_output, "Signal", "Trait_raising", "Other_allele", "Weight", "SE_weight", "a1", "a2", "beta", "se" )


# renaming the columns for ease of use 
allele_matching <- allele_matching %>%
  rename(
    SNP_T = "Signal",
    Effect_allele_T = "Trait_raising",
    Reference_allele_T = "Other_allele",
    BETA_T = "Weight",
    SE_T = "SE_weight",
    Reference_allele_FG = "a2",
    Effect_allele_FG = "a1",
    male_beta_FG = "beta",
    male_se_FG = "se"
  )


# identify trait increasing allele for T

allele_matching$T_inc_allele <- if_else(allele_matching$BETA_T<0, allele_matching$Reference_allele_T, 
                                           allele_matching$Effect_allele_T)

allele_matching$BETA_T <- as.numeric(allele_matching$BETA_T)
allele_matching$ABS_BETA_T <- abs(allele_matching$BETA_T)

# harmonising so the effect alleles for CAD and SHBG are the same
# changing the betas here 

allele_matching$male_beta_FG <- as.numeric(allele_matching$male_beta_FG)

allele_matching$HARM_MALE_BETA_FG <- if_else(allele_matching$T_inc_allele!=allele_matching$Effect_allele_FG,
                                              allele_matching$male_beta_FG*-1, allele_matching$male_beta_FG)




plot(allele_matching$ABS_BETA_T, allele_matching$HARM_MALE_BETA_FG)
allele_matching$male_se_FG <- as.numeric(allele_matching$male_se_FG)
IVW_weights <- allele_matching$male_se_FG^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_FG ~ allele_matching$ABS_BETA_T- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red")
summary_model <- summary(inverse_weighted_LR)
summary_model


MRObject = mr_input(bx = allele_matching$ABS_BETA_T, bxse = allele_matching$SE_T, 
                    by = allele_matching$HARM_MALE_BETA_FG, byse = allele_matching$male_se_FG, snps = allele_matching$SNP_T)

mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)

mr_allmethods(MRObject)
mr_plot(mr_allmethods(MRObject))
