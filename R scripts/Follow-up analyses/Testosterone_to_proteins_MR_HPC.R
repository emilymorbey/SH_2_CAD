

.libPaths("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/R_lib/")
sink("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Proteomics/Proteomics.log")

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)


#Load required packages
install.packages("tidyverse", dependencies = TRUE)
library(tidyverse) 
install.packages("data.table")
library(data.table)


setwd("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Proteomics/OL_Biobank/Males") 

files <- read.table(file="/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Proteomics/OL_Biobank/Males/Males_OLINK.csv", header=TRUE)

SNPs <- read.table(file="/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Proteomics/Testosterone_SNPs_OL.csv", header = TRUE)

# Initialize an empty dataframe to store the results
results_df <- data.frame(File = character(), Estimate = numeric(), SE = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)

# Loop over each file in the list
for(i in seq_along(files[,1])) {
  # Read the file
  T <- read.table(files[i,1], header = TRUE)
  
  # Merge with the SNP data
  merged_data <- merge(SNPs, T, by = "ID", all.x = TRUE)
  
  # Select relevant columns 
  allele_matching <- select(merged_data, "ID", "Trait_raising", 
                            "Other_allele", "Weight", "SE_weight", 
                            "MarkerName", "ALLELE0", "ALLELE1",
                            "BETA", "SE")

  allele_matching <- allele_matching[!allele_matching$ID == "rs56196860", ]
  
  # Harmonizing so the effect alleles for testosterone and protein are the same
  # Changing the betas here 


allele_matching$T_inc_allele <- if_else(allele_matching$Weight<0, allele_matching$Other_allele, 
                                           allele_matching$Trait_raising)

allele_matching$Weight <- as.numeric(allele_matching$Weight)
allele_matching$ABS_BETA_T <- abs(allele_matching$Weight)

# harmonising so the effect alleles for CAD and SHBG are the same
# changing the betas here 

allele_matching$BETA <- as.numeric(allele_matching$BETA)

allele_matching$HARM_MALE_BETA_PROTEIN <- if_else(allele_matching$T_inc_allele!=allele_matching$ALLELE1,
                                              allele_matching$BETA*-1, allele_matching$BETA)
  
  # Running the model 
  allele_matching$SE <- as.numeric(allele_matching$SE)
  IVW_weights <- allele_matching$SE^-2 
  inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_PROTEIN ~ allele_matching$ABS_BETA_T - 1, weights = IVW_weights) 
  
  # Get the summary of the model
  summary_model <- summary(inverse_weighted_LR)
  
  # Extract the estimate and p-value
  estimate <- summary_model$coefficients[1, "Estimate"]
  se <- summary_model$coefficients[1, "Std. Error"]
  p_value <- summary_model$coefficients[1, "Pr(>|t|)"]
  
  # Add the results to the dataframe
  results_df <- rbind(results_df, data.frame(File = files[i, 1], 
  Estimate = estimate, SE = se, P_Value = p_value, stringsAsFactors = FALSE))
}

# Print the results
print(results_df)

olink_info <- read.table(file="/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Proteomics/OL_Biobank/OLINK_PROTEINS_AND_IDS.csv", header = TRUE)

results_df$Olink.ID <- sub(".*_(.*?)\\..*", "\\1", results_df$File)

olink_hits <- merge(results_df, olink_info, by = "Olink.ID", all = TRUE)

olink_hits <- na.omit(olink_hits[!is.na(olink_hits$File), ])

olink_bonferroni_sig_hits <- olink_hits %>% filter(P_Value < 0.0000346)


write.table(olink_bonferroni_sig_hits, file="/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Proteomics/OL_Biobank/Males/Male_OL_Biobank_T_MR_9.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
