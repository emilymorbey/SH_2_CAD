library(tidyverse)
library(readxl)
library(MendelianRandomization)



CAD_to_T <- read_excel("not found inputs/CADtoT.xlsx", sheet = "FOUND_ON_RSID")
CAD_to_T_2 <- read_excel("not found inputs/CADtoT.xlsx", sheet = "GWAS_ESSENTIALS")

names(CAD_to_T_2)[names(CAD_to_T_2) == "rsID"] <- "SNP"

dat <- merge(CAD_to_T_2, CAD_to_T, by = "SNP" )



# renaming the columns for ease of use 
allele_matching <- dat %>%
  rename(
    SNP = "SNP",
    Effect_allele_CAD = "Effect_allele (EA)",
    Reference_allele_CAD = "Non_effect_allele",
    EA_FREQ_CAD = "EA_freq",
    BETA_CAD = "Beta",
    SE_CAD = "SE.x",
    ODDS_RATIO_CAD = "Odds_ratio",
    Effect_allele_T = "ALLELE1",
    Reference_allele_T = "ALLELE0",
    EA_FREQ_T = "A1FREQ",
    BETA_T = "BETA",
    SE_T = "SE.y", 
  )


# CHOOSING PREFERRED COLUMNS


allele_matching <- select(allele_matching, "SNP", "Effect_allele_CAD", 
                          "Reference_allele_CAD", "EA_FREQ_CAD", "BETA_CAD", 
                          "SE_CAD", "ODDS_RATIO_CAD", "Effect_allele_T",
                          "Reference_allele_T", "EA_FREQ_T", "BETA_T", 
                          "SE_T")



# identify trait increasing allele for SHBG

allele_matching$CAD_inc_allele <- if_else(allele_matching$BETA_CAD<0, allele_matching$Reference_allele_CAD, 
                                        allele_matching$Effect_allele_CAD)

allele_matching$BETA_CAD <- as.numeric(allele_matching$BETA_CAD)
allele_matching$ABS_BETA_CAD <- abs(allele_matching$BETA_CAD)

# harmonising so the effect alleles for CAD and SHBG are the same
# changing the betas here 

allele_matching$BETA_T <- as.numeric(allele_matching$BETA_T)

allele_matching$HARM_BETA_T <- if_else(allele_matching$CAD_inc_allele!=allele_matching$Effect_allele_T,
                                              allele_matching$BETA_T*-1, allele_matching$BETA_T)



# PLOTTING IT 

plot(allele_matching$ABS_BETA_CAD, allele_matching$HARM_BETA_T)
allele_matching$SE_T <- as.numeric(allele_matching$SE_T)
IVW_weights <- allele_matching$SE_T^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_BETA_T ~ allele_matching$ABS_BETA_CAD- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red")
summary_model <- summary(inverse_weighted_LR)
summary_model





MRObject = mr_input(bx = allele_matching$ABS_BETA_CAD, bxse = allele_matching$SE_CAD, 
                    by = allele_matching$HARM_BETA_T, byse = allele_matching$SE_T, snps = allele_matching$SNP)


mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)
