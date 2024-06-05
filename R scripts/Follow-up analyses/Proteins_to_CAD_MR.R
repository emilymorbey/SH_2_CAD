library(tidyverse)
library(readxl)

setwd("C:/Users/emorb/OneDrive - University of Cambridge/PhD/MR/Testosterone_CAD_MR/Testosterone CAD MR R files/TestosteroneCAD/Proteomics")
allele_matching <- read_excel("PROTEIN_TOP_HITS_AND_CAD.xlsx", sheet = "T&P E&O")




allele_matching <- select(allele_matching, "OLID", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "SE", "reference_allele", "other_allele", "eaf", "male_beta", "male_se" )

# renaming the columns for ease of use 
allele_matching <- allele_matching %>%
  rename(
    OLID = "OLID",
    Effect_allele_P = "ALLELE1",
    Reference_allele_P = "ALLELE0",
    EA_FREQ_P = "A1FREQ",
    BETA_P = "BETA",
    SE_P = "SE",
    Reference_allele_CAD = "reference_allele",
    Effect_allele_CAD = "other_allele",
    EA_FREQ_CAD = "eaf",
    male_beta_CAD = "male_beta",
    male_se_CAD = "male_se"
  )


# identify trait increasing allele for SHBG

allele_matching$P_inc_allele <- if_else(allele_matching$BETA_P<0, allele_matching$Reference_allele_P, 
                                        allele_matching$Effect_allele_P)

allele_matching$BETA_P <- as.numeric(allele_matching$BETA_P)
allele_matching$ABS_BETA_P <- abs(allele_matching$BETA_P)

# harmonising so the effect alleles for CAD and SHBG are the same
# changing the betas here 

allele_matching$male_beta_CAD <- as.numeric(allele_matching$male_beta_CAD)

allele_matching$HARM_MALE_BETA_CAD <- if_else(allele_matching$P_inc_allele!=allele_matching$Effect_allele_CAD,
                                              allele_matching$male_beta_CAD*-1, allele_matching$male_beta_CAD)


allele_matching$male_se_CAD <- as.numeric(allele_matching$male_se_CAD)


OL20090 <- allele_matching %>% 
  filter(OLID=="OL20090")

OL20204 <- allele_matching %>% 
  filter(OLID=="OL20204")

OL20220 <- allele_matching %>% 
  filter(OLID=="OL20220")

OL20245 <- allele_matching %>% 
  filter(OLID=="OL20245")

OL20758 <- allele_matching %>% 
  filter(OLID=="OL20758")

OL21033 <- allele_matching %>% 
  filter(OLID=="OL21033")

OL21257 <- allele_matching %>% 
  filter(OLID=="OL21257")

OL20967 <- allele_matching %>% 
  filter(OLID=="OL20967")


 
 library(dplyr)
 

 # Function to compute IVW estimate, standard error, OR, LCI, and UCI
 compute_ivw <- function(data, olid) {
   beta.ivw <- sum(data$BETA_P * data$male_beta_CAD * data$male_se_CAD^-2) / sum(data$BETA_P^2 * data$male_se_CAD^-2)
   se.ivw <- 1 / sqrt(sum(data$BETA_P^2 * data$male_se_CAD^-2))
   or <- exp(beta.ivw)
   lci <- exp(beta.ivw - (2 * se.ivw))
   uci <- exp(beta.ivw + (2 * se.ivw))
   return(data.frame(OLID = olid, beta = beta.ivw, se = se.ivw, OR = or, LCI = lci, UCI = uci))
 }
 
 # List of OLIDs to filter and process
 olid_list <- c("OL20090", "OL20204", "OL20220", "OL20245", "OL20758", "OL21033", "OL21257", "OL20967")
 
 # Initialize an empty dataframe to store the results
 results <- data.frame(OLID = character(), beta = numeric(), se = numeric(), OR = numeric(), LCI = numeric(), UCI = numeric(), stringsAsFactors = FALSE)
 
 # Iterate over each OLID, filter the data, and compute IVW estimate
 for (olid in olid_list) {
   df <- allele_matching %>% filter(OLID == olid)
   result <- compute_ivw(df, olid)
   results <- rbind(results, result)
 }
 
 # Print the results
 print(results)
 
 
write.csv(results, "ivw_results_protein_cad.csv", row.names = FALSE)