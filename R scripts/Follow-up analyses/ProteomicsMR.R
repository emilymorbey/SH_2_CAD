library(tidyverse)
library(MendelianRandomization)
library(data.table)
library(readr)


setwd("C:/Users/emorb/OneDrive - University of Cambridge/PhD/MR/Testosterone_CAD_MR/Testosterone CAD MR R files/TestosteroneCAD/Proteomics")


files <- read.table("Males_OLINK.csv", header=TRUE)

SNPs <- read.csv("Testosterone_SNPs_OL.csv", header = TRUE)

# Initialize an empty dataframe to store the results
results_df <- data.frame(File = character(), Estimate = numeric(), SE = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)

# Loop over each file in the list
for(i in seq_along(files[,1])) {
  # Read the file
  T <- read.table(files[i,1], header = TRUE)
  
  # Merge with the SNP data
  merged_data <- merge(SNPs, T, by = "ID")
  
  # Select relevant columns 
  allele_matching <- select(merged_data, "ID", "Trait_raising", 
                            "Other_allele", "Weight", "SE_weight", 
                            "MarkerName", "ALLELE0", "ALLELE1",
                            "BETA", "SE")
  
  # Harmonizing so the effect alleles for CAD and SHBG are the same
  # Changing the betas here 
  allele_matching$BETA <- as.numeric(allele_matching$BETA)
  allele_matching$HARM_BETA <- if_else(allele_matching$Trait_raising != allele_matching$ALLELE0,
                                       allele_matching$BETA * -1, allele_matching$BETA)
  
  # Running the model 
  allele_matching$SE <- as.numeric(allele_matching$SE)
  IVW_weights <- allele_matching$SE^-2 
  inverse_weighted_LR <- lm(allele_matching$HARM_BETA ~ allele_matching$Weight - 1, weights = IVW_weights) 
  
  # Get the summary of the model
  summary_model <- summary(inverse_weighted_LR)
  
  # Extract the estimate and p-value
  estimate <- summary_model$coefficients[1, "Estimate"]
  se <- summary_model$coefficients[1, "Std. Error"]
  p_value <- summary_model$coefficients[1, "Pr(>|t|)"]
  
  # Add the results to the dataframe
  results_df <- rbind(results_df, data.frame(File = files[i, 1], Estimate = estimate, SE = se, P_Value = p_value, stringsAsFactors = FALSE))
}

# Print the results
print(results_df)




olinkid <- read.csv("OLINKIDS.csv", header = TRUE)
olinkname <- read.csv("OLINKNAMES.csv", header = TRUE)


olink <- merge(olinkid, olinkname, by = "ID")

files <- read.delim("Male_OL_Biobank_T_MR.txt")


files$Olink.ID <- sub(".*_(.*?)\\..*", "\\1", files$File)

olink_hits <- merge(files, olink, by = "Olink.ID")

olink_med_hits <- olink_hits %>% filter(P_Value < 0.00005)
olink_sig_hits <- olink_hits %>% filter(P_Value < 0.00000005)



write.table(olink_med_hits, "olinkmedhits.tsv", sep = "\t", row.names = TRUE)
write.table(olink_sig_hits, "olinksighits.tsv", sep = "\t", row.names = TRUE)        
