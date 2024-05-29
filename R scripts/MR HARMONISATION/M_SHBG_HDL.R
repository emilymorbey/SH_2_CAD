library(tidyverse)
library(readxl)
library(MendelianRandomization)

##################################################################################

# HARMONISATION AND MR

##################################################################################

# looking at the allele matching and frequencies etc.
M_SHBG_proxies_output <- read_excel("not found inputs/SNPs_M_SHBG_AND_HDL.xlsx", sheet = "T&P E&O")



allele_matching <- select(M_SHBG_proxies_output, "SNP.x", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE.x", "REF", "ALT", "POOLED_ALT_AF", "EFFECT_SIZE", "SE.y" )

# renaming the columns for ease of use 
allele_matching <- allele_matching %>%
  rename(
    SNP_SHBG = "SNP.x",
    other_allele_SHBG = "ALLELE1",
    reference_allele_SHBG = "ALLELE0",
    A1FREQ_SHBG = "A1FREQ",
    BETA_SHBG = "BETA",
    SE_SHBG = "SE.x",
    reference_allele_HDL = "REF",
    other_allele_HDL = "ALT",
    eaf_HDL = "POOLED_ALT_AF",
    beta_HDL = "EFFECT_SIZE",
    se_HDL = "SE.y"
  )

# identify trait increasing allele for SHBG

allele_matching$SHBG_inc_allele <- if_else(allele_matching$BETA_SHBG<0, allele_matching$reference_allele_SHBG, 
                                           allele_matching$other_allele_SHBG)

allele_matching$BETA_SHBG <- as.numeric(allele_matching$BETA_SHBG)
allele_matching$ABS_BETA_SHBG <- abs(allele_matching$BETA_SHBG)

# harmonising so the effect alleles for CAD and SHBG are the same
# changing the betas here 

allele_matching$beta_HDL <- as.numeric(allele_matching$beta_HDL)

allele_matching$HARM_BETA_HDL <- if_else(allele_matching$SHBG_inc_allele!=allele_matching$other_allele_HDL,
                                         allele_matching$beta_HDL*-1, allele_matching$beta_HDL)




plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_BETA_HDL)


M_SHBG_proxies_output$SE.y <- as.numeric(M_SHBG_proxies_output$SE.y)
IVW_weights <- M_SHBG_proxies_output$SE.y^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_BETA_HDL ~ allele_matching$ABS_BETA_SHBG- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red")
summary_model <- summary(inverse_weighted_LR)
summary_model



allele_matching$ABS_BETA_SHBG <- as.numeric(allele_matching$ABS_BETA_SHBG)
allele_matching$se_HDL <- as.numeric(allele_matching$se_HDL)
allele_matching$SE_SHBG <- as.numeric(allele_matching$SE_SHBG)



MRObject = mr_input(bx = allele_matching$ABS_BETA_SHBG, bxse = allele_matching$SE_SHBG, 
                    by = allele_matching$HARM_BETA_HDL, byse = allele_matching$se_HDL, snps = allele_matching$SNP_SHBG)

mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)


mr_allmethods(MRObject)
mr_plot(mr_allmethods(MRObject))



plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_BETA_HDL  ,
     xlab = "SNP effect on SHBG",  # Replace with your desired x-axis label
     ylab = "SNP effect on HDL-c",
     main = "Male SHBG")  # Replace with your desired y-axis label

M_SHBG_proxies_output$SE.y <- as.numeric(M_SHBG_proxies_output$SE.y)
IVW_weights <- M_SHBG_proxies_output$SE.y^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_BETA_HDL ~ allele_matching$ABS_BETA_SHBG - 1, weights = IVW_weights)
summary(inverse_weighted_LR)

abline(inverse_weighted_LR, col = "red")




SHBG_SNPS <- allele_matching %>%
  filter(SNP_SHBG=="rs6258")

snp_x <- 0.610432
snp_y <- 0.0196715
snp_label <- "rs6258"

plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_BETA_HDL,
     xlab = "SNP effect on SHBG",
     ylab = "SNP effect on HDL-c",
     main = "Male SHBG")

M_SHBG_proxies_output$SE.y <- as.numeric(M_SHBG_proxies_output$SE.y)
IVW_weights <- M_SHBG_proxies_output$SE.y^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_BETA_HDL ~ allele_matching$ABS_BETA_SHBG - 1, weights = IVW_weights)
summary(inverse_weighted_LR)

abline(inverse_weighted_LR, col = "red")


# Additional code to label the SNP
text(snp_x, snp_y, snp_label, col = "blue", pos = 1, cex = 0.7)


#######################################################################################

### other interprative graphs #############

mr_plot(MRObject, interactive=FALSE, labels=TRUE)
mr_forest(MRObject, ordered=TRUE)
mr_loo(MRObject)
mr_funnel(MRObject)
