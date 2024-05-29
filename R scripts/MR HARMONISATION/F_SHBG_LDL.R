library(tidyverse)
library(readxl)
library(MendelianRandomization)

##################################################################################

# HARMONISATION AND MR

##################################################################################

# looking at the allele matching and frequencies etc.
F_SHBG_proxies_output <- read_excel("not found inputs/SNPs_F_SHBG_AND_LDL.xlsx", sheet = "T&P E&O")



allele_matching <- select(F_SHBG_proxies_output, "SNP.x", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE.x", "REF", "ALT", "POOLED_ALT_AF", "EFFECT_SIZE", "SE.y" )

# renaming the columns for ease of use 
allele_matching <- allele_matching %>%
  rename(
    SNP_SHBG = "SNP.x",
    other_allele_SHBG = "ALLELE1",
    reference_allele_SHBG = "ALLELE0",
    A1FREQ_SHBG = "A1FREQ",
    BETA_SHBG = "BETA",
    SE_SHBG = "SE.x",
    reference_allele_LDL = "REF",
    other_allele_LDL = "ALT",
    eaf_LDL = "POOLED_ALT_AF",
    beta_LDL = "EFFECT_SIZE",
    se_LDL = "SE.y"
  )

# identify trait increasing allele for SHBG

allele_matching$SHBG_inc_allele <- if_else(allele_matching$BETA_SHBG<0, allele_matching$reference_allele_SHBG, 
                                        allele_matching$other_allele_SHBG)

allele_matching$BETA_SHBG <- as.numeric(allele_matching$BETA_SHBG)
allele_matching$ABS_BETA_SHBG <- abs(allele_matching$BETA_SHBG)

# harmonising so the effect alleles for CAD and SHBG are the same
# changing the betas here 

allele_matching$beta_LDL <- as.numeric(allele_matching$beta_LDL)

allele_matching$HARM_BETA_LDL <- if_else(allele_matching$SHBG_inc_allele!=allele_matching$other_allele_LDL,
                                         allele_matching$beta_LDL*-1, allele_matching$beta_LDL)




plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_BETA_LDL)


F_SHBG_proxies_output$SE.y <- as.numeric(F_SHBG_proxies_output$SE.y)
IVW_weights <- F_SHBG_proxies_output$SE.y^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_BETA_LDL ~ allele_matching$ABS_BETA_SHBG- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red")
summary_model <- summary(inverse_weighted_LR)
summary_model



allele_matching$ABS_BETA_SHBG <- as.numeric(allele_matching$ABS_BETA_SHBG)
allele_matching$se_LDL <- as.numeric(allele_matching$se_LDL)
allele_matching$SE_SHBG <- as.numeric(allele_matching$SE_SHBG)



MRObject = mr_input(bx = allele_matching$ABS_BETA_SHBG, bxse = allele_matching$SE_SHBG, 
                    by = allele_matching$HARM_BETA_LDL, byse = allele_matching$se_LDL, snps = allele_matching$SNP_SHBG)

mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)





plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_BETA_LDL  ,
     xlab = "SNP effect on SHBG",  # Replace with your desired x-axis label
     ylab = "SNP effect on LDL-c",
     main = "Female SHBG")  # Replace with your desired y-axis label

F_SHBG_proxies_output$SE.y <- as.numeric(F_SHBG_proxies_output$SE.y)
IVW_weights <- F_SHBG_proxies_output$SE.y^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_BETA_LDL ~ allele_matching$ABS_BETA_SHBG - 1, weights = IVW_weights)
summary(inverse_weighted_LR)

abline(inverse_weighted_LR, col = "red")



SHBG_SNPS <- allele_matching %>%
  filter(SNP_SHBG=="rs6258")

snp_x <- 0.67
snp_y <- -0.044
snp_label <- "rs6258"

plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_BETA_LDL,
     xlab = "SNP effect on SHBG",
     ylab = "SNP effect on LDL-c",
     main = "Female SHBG")

F_SHBG_proxies_output$SE.y <- as.numeric(F_SHBG_proxies_output$SE.y)
IVW_weights <- F_SHBG_proxies_output$SE.y^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_BETA_LDL ~ allele_matching$ABS_BETA_SHBG - 1, weights = IVW_weights)
summary(inverse_weighted_LR)

abline(inverse_weighted_LR, col = "red")


# Additional code to label the SNP
text(snp_x, snp_y, snp_label, col = "blue", pos = 1, cex = 0.7)

##################################################################################




#######################################################################################

### other interprative graphs #############

mr_plot(MRObject, interactive=FALSE, labels=TRUE)
mr_forest(MRObject, ordered=TRUE)
mr_loo(MRObject)
mr_funnel(MRObject)
