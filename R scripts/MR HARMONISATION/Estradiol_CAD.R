library(tidyverse)
library(readxl)
library(MendelianRandomization)

#################################################################################
# don't run this again - it was run once and that was fine
# scroll down

# first read in the original GWAS output for the hormone
# second in the file which has the SNPs that were found in the CAD file from the
# hormone file


M_E <- read_excel("Not in CAD inputs/Original Files/M_estradiol.xlsx")
M_E_AND_CAD <- read_excel("Not in CAD inputs/SNPs in both datasets/SNPs_M_Estradiol_AND_CAD.xlsx", sheet = "RSID+CHRPOS")


M_E <- M_E %>% 
  rename(
    rsid_ukb = Signal,
  )

merged_df <- merge(M_E, M_E_AND_CAD, by.x = 'rsid_ukb', by.y = 'rsid_ukb', all = TRUE)

sum(is.na(merged_df$rs_number))

not_in_CAD <- subset(merged_df, is.na(rs_number))

# writing this dataframe into an excel spreadsheet 

library(openxlsx)

write.xlsx(not_in_CAD, "not_in_cad_M_Estradiol.xlsx", rowNames = TRUE)







##################################################################################

# HARMONISATION AND MR

##################################################################################

# looking at the allele matching and frequencies etc.
M_E_proxies_output <- read_excel("not found inputs/SNPs_M_Estradiol_AND_CAD.xlsx", sheet = "T&P E&O")



allele_matching <- select(M_E_proxies_output, "SNP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "reference_allele", "other_allele", "eaf", "male_beta", "male_se" )

# renaming the columns for ease of use 
allele_matching <- allele_matching %>%
  rename(
    SNP_E = "SNP",
    ALLELE1_E = "ALLELE1",
    ALLELE0_E = "ALLELE0",
    A1FREQ_E = "A1FREQ",
    BETA_E = "BETA",
    SE_E = "SE",
    reference_allele_CAD = "reference_allele",
    other_allele_CAD = "other_allele",
    eaf_CAD = "eaf",
    male_beta_CAD = "male_beta",
    male_se_CAD = "male_se"
  )

# identify trait increasing allele for SHBG

allele_matching$E_inc_allele <- if_else(allele_matching$BETA_E<0, allele_matching$ALLELE0_E, 
                                        allele_matching$ALLELE1_E)

allele_matching$BETA_E <- as.numeric(allele_matching$BETA_E)
allele_matching$ABS_BETA_E <- abs(allele_matching$BETA_E)

# harmonising so the effect alleles for CAD and SHBG are the same
# changing the betas here 

allele_matching$male_beta_CAD <- as.numeric(allele_matching$male_beta_CAD)

allele_matching$HARM_MALE_BETA_CAD <- if_else(allele_matching$E_inc_allele!=allele_matching$other_allele_CAD,
                                                allele_matching$male_beta_CAD*-1, allele_matching$male_beta_CAD)




plot(allele_matching$ABS_BETA_E, allele_matching$HARM_MALE_BETA_CAD)
M_E_proxies_output$male_se <- as.numeric(M_E_proxies_output$male_se)
IVW_weights <- M_E_proxies_output$male_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_E- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red")
summary_model <- summary(inverse_weighted_LR)
summary_model



M_E_proxies_output$male_beta <- as.numeric(M_E_proxies_output$male_beta)

allele_matching$ABS_BETA_E <- as.numeric(allele_matching$ABS_BETA_E)
allele_matching$male_se_CAD <- as.numeric(allele_matching$male_se_CAD)
allele_matching$SE_E <- as.numeric(allele_matching$SE_E)

MRObject = mr_input(bx = allele_matching$ABS_BETA_E, bxse = allele_matching$SE_E, 
                    by = allele_matching$HARM_MALE_BETA_CAD, byse = allele_matching$male_se_CAD, snps = allele_matching$SNP_E)

mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)


mr_allmethods(MRObject)
mr_plot(mr_allmethods(MRObject))



plot(allele_matching$ABS_BETA_E, allele_matching$HARM_MALE_BETA_CAD,
     xlab = "SNP effect on Estradiol",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Male Estradiol")  # Replace with your desired y-axis label

M_E_proxies_output$male_se <- as.numeric(M_E_proxies_output$male_se)
IVW_weights <- M_E_proxies_output$male_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_E - 1, weights = IVW_weights)
summary(inverse_weighted_LR)

abline(inverse_weighted_LR, col = "red")





##################################################################################




mr_plot(MRObject, interactive=FALSE, labels=TRUE)
mr_forest(MRObject, ordered=TRUE)
mr_loo(MRObject)
mr_funnel(MRObject)






















