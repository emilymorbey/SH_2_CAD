library(devtools)


devtools::install_github("rondolab/MR-PRESSO", force = TRUE)
library(MRPRESSO)

# run the M_TESTSOTERONE_CAD script before running this 

allele_matching <- as.data.frame(allele_matching)

mr_presso(BetaOutcome = "male_beta_CAD", BetaExposure = "BETA_T", SdOutcome = "male_se_CAD", SdExposure = "SE_T", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = allele_matching, NbDistribution = 3000,  SignifThreshold = 0.05)


