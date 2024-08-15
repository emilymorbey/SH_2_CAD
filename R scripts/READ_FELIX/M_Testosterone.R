##################################################################################

# HARMONISATION AND MR

##################################################################################
library(tidyverse)
library(readxl)
library(MendelianRandomization)


setwd("C:/Users/emorb/OneDrive - University of Cambridge/PhD/MR/Testosterone_CAD_MR/Testosterone CAD MR R files")

M_T_proxies_output <- read_excel("TestosteroneCAD/not found inputs/SNPs_M_Testosterone_AND_CAD.xlsx", sheet = "FREE T E&O")

### selecting appropriate columns for harmonisation 

allele_matching <- select(M_T_proxies_output, "SNP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "reference_allele", "other_allele", "eaf", "male_beta", "male_se" )

# renaming the columns for ease of use 
allele_matching <- allele_matching %>%
  rename(
    SNP_T = "SNP",
    Effect_allele_T = "ALLELE1",
    Reference_allele_T = "ALLELE0",
    EA_FREQ_T = "A1FREQ",
    BETA_T = "BETA",
    SE_T = "SE",
    Reference_allele_CAD = "reference_allele",
    Effect_allele_CAD = "other_allele",
    EA_FREQ_CAD = "eaf",
    male_beta_CAD = "male_beta",
    male_se_CAD = "male_se"
  )

# identify trait increasing allele for Testosterone
# here we are saying, if beta is negative, then the reference allele is the trait increasing allele, if beta is positive, then the effect allele is the trait increasing allele

allele_matching$T_inc_allele <- if_else(allele_matching$BETA_T<0, allele_matching$Reference_allele_T, 
                                           allele_matching$Effect_allele_T)



### setting the betas as numeric 

allele_matching$BETA_T <- as.numeric(allele_matching$BETA_T)
allele_matching$ABS_BETA_T <- abs(allele_matching$BETA_T)

# harmonising so the effect alleles for CAD and Testosterone are the same
# changing the betas of the CAD SNPs to match the new effect allele for CAD

allele_matching$male_beta_CAD <- as.numeric(allele_matching$male_beta_CAD)

# here we are saying, if the trait increasing allele for Testosterone is not the same as the effect allele for CAD, then the beta for CAD is multiplied by -1, otherwise it is the same

allele_matching$HARM_MALE_BETA_CAD <- if_else(allele_matching$T_inc_allele!=allele_matching$Effect_allele_CAD,
                                              allele_matching$male_beta_CAD*-1, allele_matching$male_beta_CAD)



## manually running the IVW method

plot(allele_matching$ABS_BETA_T, allele_matching$HARM_MALE_BETA_CAD)
allele_matching$male_se_CAD <- as.numeric(allele_matching$male_se_CAD)
IVW_weights <- allele_matching$male_se_CAD^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_T- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red")
summary_model <- summary(inverse_weighted_LR)
summary_model

### making some more things numeric 

M_T_proxies_output$male_beta <- as.numeric(M_T_proxies_output$male_beta)
allele_matching$ABS_BETA_T <- as.numeric(allele_matching$ABS_BETA_T)
allele_matching$male_se_CAD <- as.numeric(allele_matching$male_se_CAD)
allele_matching$SE_T <- as.numeric(allele_matching$SE_T)

### running all models using the MR package

MRObject = mr_input(bx = allele_matching$ABS_BETA_T, bxse = allele_matching$SE_T, 
                    by = allele_matching$HARM_MALE_BETA_CAD, byse = allele_matching$male_se_CAD, snps = allele_matching$SNP_T)

mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)

mr_allmethods(MRObject)
mr_plot(mr_allmethods(MRObject))


### creating a plot for the MR object with a line for the IVW method

plot(allele_matching$ABS_BETA_T, allele_matching$HARM_MALE_BETA_CAD,
     xlab = "SNP effect on Testosterone",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Male Testosterone")  # Replace with your desired y-axis label

M_T_proxies_output$male_se <- as.numeric(M_T_proxies_output$male_se)
IVW_weights <- M_T_proxies_output$male_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_T - 1, weights = IVW_weights)
summary(inverse_weighted_LR)

abline(inverse_weighted_LR, col = "red")





#######################################################################################

### other interprative graphs #############

mr_plot(MRObject, interactive=TRUE, labels=TRUE)
mr_forest(MRObject, ordered=TRUE)
mr_loo(MRObject)
plot_object <- mr_loo(MRObject)


customized_plot <- plot_object +
  ggtitle("Leave-One-Out MR Plot") +  # Adding a title
  xlab("Effect size estimate") +                       # Custom x-axis label
  ylab("SNPs") +               # Custom y-axis label
  theme_minimal() +                   # Applying a minimal theme for a clean look
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5), # Centered title
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",              # Moving the legend to the bottom
    legend.title = element_blank(),          # Removing the legend title for simplicity
    axis.text.y = element_blank(),           # Removing y-axis text
    axis.title.y = element_blank(),          # Removing y-axis title
    axis.ticks.y = element_blank(),          # Removing y-axis ticks
    panel.grid.major = element_blank(),  # Light grey grid lines
    panel.grid.minor = element_blank(),                  # Removing minor grid lines
    plot.background = element_rect(fill = "white", color = "white"), # White background
    panel.border = element_blank()           # Removing panel border
  ) +
  geom_point(color = "blue", linewidth = 2) +        # Blue points for better visibility
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1) # Red trend line

# Display the customized plot
print(customized_plot)

mr_funnel(MRObject)
??mr_funnel


mr_funnel(MRObject) + ggtitle("Funnel plot for male testosterone and CAD")


################ RUNNING WITHOUT OUTLIER ############################



# leaving out the SNP that is causing the problem

allele_matching <- allele_matching[!allele_matching$SNP_T == "rs56196860", ]
M_T_proxies_output <- M_T_proxies_output[!M_T_proxies_output$Target == "rs56196860", ]
# identify trait increasing allele for SHBG

allele_matching$T_inc_allele <- if_else(allele_matching$BETA_T<0, allele_matching$Reference_allele_T, 
                                        allele_matching$Effect_allele_T)

allele_matching$BETA_T <- as.numeric(allele_matching$BETA_T)
allele_matching$ABS_BETA_T <- abs(allele_matching$BETA_T)

# harmonising so the effect alleles for CAD and SHBG are the same
# changing the betas here 

allele_matching$male_beta_CAD <- as.numeric(allele_matching$male_beta_CAD)

allele_matching$HARM_MALE_BETA_CAD <- if_else(allele_matching$T_inc_allele!=allele_matching$Effect_allele_CAD,
                                              allele_matching$male_beta_CAD*-1, allele_matching$male_beta_CAD)




plot(allele_matching$ABS_BETA_T, allele_matching$HARM_MALE_BETA_CAD)
M_T_proxies_output$male_se <- as.numeric(M_T_proxies_output$male_se)
IVW_weights <- M_T_proxies_output$male_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_T- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red")
summary_model <- summary(inverse_weighted_LR)
summary_model



M_T_proxies_output$male_beta <- as.numeric(M_T_proxies_output$male_beta)

allele_matching$ABS_BETA_T <- as.numeric(allele_matching$ABS_BETA_T)
allele_matching$male_se_CAD <- as.numeric(allele_matching$male_se_CAD)
allele_matching$SE_T <- as.numeric(allele_matching$SE_T)

MRObject = mr_input(bx = allele_matching$ABS_BETA_T, bxse = allele_matching$SE_T, 
                    by = allele_matching$HARM_MALE_BETA_CAD, byse = allele_matching$male_se_CAD, snps = allele_matching$SNP_T)

mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)

mr_allmethods(MRObject)
mr_plot(mr_allmethods(MRObject))
mr_plot(MRObject, interactive=TRUE, labels=TRUE)


plot(allele_matching$ABS_BETA_T, allele_matching$HARM_MALE_BETA_CAD,
     xlab = "SNP effect on Testosterone",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Male Testosterone")  # Replace with your desired y-axis label

M_T_proxies_output$male_se <- as.numeric(M_T_proxies_output$male_se)
IVW_weights <- M_T_proxies_output$male_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_T - 1, weights = IVW_weights)
summary(inverse_weighted_LR)

abline(inverse_weighted_LR, col = "red")




#### USING WEIGHTS FROM CLUSTER 


##################################################################################

# HARMONISATION AND MR

##################################################################################
library(tidyverse)
library(readxl)
library(MendelianRandomization)


setwd("C:/Users/emorb/OneDrive - University of Cambridge/PhD/MR/Testosterone_CAD_MR/Testosterone CAD MR R files")

M_T_proxies_output <- read_excel("TestosteroneCAD/not found inputs/SNPs_M_Testosterone_AND_CAD.xlsx", sheet = "all 3")
M_T_proxies_output <- M_T_proxies_output[!M_T_proxies_output$Signal == "rs56196860", ]

head(M_T_proxies_output)
View(M_T_proxies_output)
### selecting appropriate columns for harmonisation 

allele_matching <- select(M_T_proxies_output, "Signal", "Trait_raising", "Other_allele", "Weight", "SE_weight", "reference_allele", "other_allele", "male_beta", "male_se" )

# renaming the columns for ease of use 
allele_matching <- allele_matching %>%
  rename(
    SNP_T = "Signal",
    Effect_allele_T = "Trait_raising",
    Reference_allele_T = "Other_allele",
    BETA_T = "Weight",
    SE_T = "SE_weight",
    Reference_allele_CAD = "reference_allele",
    Effect_allele_CAD = "other_allele",
    male_beta_CAD = "male_beta",
    male_se_CAD = "male_se"
  )

# identify trait increasing allele for Testosterone
# here we are saying, if beta is negative, then the reference allele is the trait increasing allele, if beta is positive, then the effect allele is the trait increasing allele

allele_matching$T_inc_allele <- if_else(allele_matching$BETA_T<0, allele_matching$Reference_allele_T, 
                                           allele_matching$Effect_allele_T)



### setting the betas as numeric 

allele_matching$BETA_T <- as.numeric(allele_matching$BETA_T)
allele_matching$ABS_BETA_T <- abs(allele_matching$BETA_T)

# harmonising so the effect alleles for CAD and Testosterone are the same
# changing the betas of the CAD SNPs to match the new effect allele for CAD

allele_matching$male_beta_CAD <- as.numeric(allele_matching$male_beta_CAD)

# here we are saying, if the trait increasing allele for Testosterone is not the same as the effect allele for CAD, then the beta for CAD is multiplied by -1, otherwise it is the same

allele_matching$HARM_MALE_BETA_CAD <- if_else(allele_matching$T_inc_allele!=allele_matching$Effect_allele_CAD,
                                              allele_matching$male_beta_CAD*-1, allele_matching$male_beta_CAD)



## manually running the IVW method

plot(allele_matching$ABS_BETA_T, allele_matching$HARM_MALE_BETA_CAD)
allele_matching$male_se_CAD <- as.numeric(allele_matching$male_se_CAD)
IVW_weights <- allele_matching$male_se_CAD^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_T- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red")
summary_model <- summary(inverse_weighted_LR)
summary_model

### making some more things numeric 

M_T_proxies_output$male_beta <- as.numeric(M_T_proxies_output$male_beta)
allele_matching$ABS_BETA_T <- as.numeric(allele_matching$ABS_BETA_T)
allele_matching$male_se_CAD <- as.numeric(allele_matching$male_se_CAD)
allele_matching$SE_T <- as.numeric(allele_matching$SE_T)

### running all models using the MR package

MRObject = mr_input(bx = allele_matching$ABS_BETA_T, bxse = allele_matching$SE_T, 
                    by = allele_matching$HARM_MALE_BETA_CAD, byse = allele_matching$male_se_CAD, snps = allele_matching$SNP_T)

mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)

mr_allmethods(MRObject)
mr_plot(mr_allmethods(MRObject))

