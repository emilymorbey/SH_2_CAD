# looking at the allele matching and frequencies etc.
M_SHBG_proxies_output <- read_excel("not found inputs/SNPs_M_SHBG_AND_CAD.xlsx", sheet = "T&P E&O")
M_SHBG_proxies_output <- M_SHBG_proxies_output[-1,]


allele_matching <- select(M_SHBG_proxies_output, "SNP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "reference_allele", "other_allele", "eaf", "male_beta", "male_se" )

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
    male_beta_CAD = "male_beta",
    male_se_CAD = "male_se"
  )

# identify trait increasing allele for SHBG

allele_matching$SHBG_inc_allele <- if_else(allele_matching$BETA_SHBG<0, allele_matching$ALLELE0_SHBG, 
                                           allele_matching$ALLELE1_SHBG)

allele_matching$BETA_SHBG <- as.numeric(allele_matching$BETA_SHBG)
allele_matching$ABS_BETA_SHBG <- abs(allele_matching$BETA_SHBG)

# harmonising so the effect alleles for CAD and SHBG are the same
# changing the betas here 

allele_matching$male_beta_CAD <- as.numeric(allele_matching$male_beta_CAD)

allele_matching$HARM_MALE_BETA_CAD <- if_else(allele_matching$SHBG_inc_allele!=allele_matching$other_allele_CAD,
                                              allele_matching$male_beta_CAD*-1, allele_matching$male_beta_CAD)



allele_matching <- filter(allele_matching, SNP_SHBG=="rs1799941")

# Scatter plot
plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_MALE_BETA_CAD, main = "SHBG gene", xlab = "SHBG", ylab = "CAD")

# Add labels to the points based on the SNP allele
for (i in seq_along(allele_matching$SNP_SHBG)) {
  text(allele_matching$ABS_BETA_SHBG[i], allele_matching$HARM_MALE_BETA_CAD[i], labels = allele_matching$SNP_SHBG[i], pos = 1)
}







plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_MALE_BETA_CAD)
M_SHBG_proxies_output$male_se <- as.numeric(M_SHBG_proxies_output$male_se)
IVW_weights <- M_SHBG_proxies_output$male_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_SHBG- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red")
summary_model <- summary(inverse_weighted_LR)
summary_model






M_SHBG_proxies_output$male_beta <- as.numeric(M_SHBG_proxies_output$male_beta)

allele_matching$ABS_BETA_SHBG <- as.numeric(allele_matching$ABS_BETA_SHBG)
allele_matching$male_se_CAD <- as.numeric(allele_matching$male_se_CAD)
allele_matching$SE_SHBG <- as.numeric(allele_matching$SE_SHBG)

MRObject = mr_input(bx = allele_matching$ABS_BETA_SHBG, bxse = allele_matching$SE_SHBG, 
                    by = allele_matching$HARM_MALE_BETA_CAD, byse = allele_matching$male_se_CAD)

mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)

mr_plot(MRObject, orientate=TRUE, line="egger")
mr_plot(mr_allmethods(MRObject, method="ivw"))


mr_allmethods(MRObject)
mr_plot(mr_allmethods(MRObject))








plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_MALE_BETA_CAD,
     xlab = "SNP effect on SHBG",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Male SHBG")  # Replace with your desired y-axis label

M_SHBG_proxies_output$male_se <- as.numeric(M_SHBG_proxies_output$male_se)
IVW_weights <- M_SHBG_proxies_output$male_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_SHBG - 1, weights = IVW_weights)
summary(inverse_weighted_LR)

abline(inverse_weighted_LR, col = "red")

