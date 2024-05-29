
### THIS IS THE CODE FOR MAKING THE SCATTER PLOTS OF THE IVW AND 
### MR EGGER ESTIMATES. YOU MUST RUN THE CODE OF THE RESPECTIVE HORMONE
### PRIOR TO RUNNING THEIR PARAGRAPH IN THIS SCRIPT




################ FEMALE SHBG ####################################

# setting up the basic plot ########################################


plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_FEMALE_BETA_CAD, pch = 16, cex = 0.7,
     xlab = "SNP effect on SHBG",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Female SHBG")

# Add error bars
segments(
  x0 = allele_matching$ABS_BETA_SHBG,
  y0 = allele_matching$HARM_FEMALE_BETA_CAD - allele_matching$female_se_CAD,
  x1 = allele_matching$ABS_BETA_SHBG,
  y1 = allele_matching$HARM_FEMALE_BETA_CAD + allele_matching$female_se_CAD,
  col = "black"
)

segments(
  x0 = allele_matching$ABS_BETA_SHBG - allele_matching$SE_SHBG, 
  y0 = allele_matching$HARM_FEMALE_BETA_CAD,
  x1 = allele_matching$ABS_BETA_SHBG + allele_matching$SE_SHBG, 
  y1 = allele_matching$HARM_FEMALE_BETA_CAD,
 col = "black"
)

# adding the lines of the different models ############################
# IVW

F_SHBG_proxies_output$female_se <- as.numeric(F_SHBG_proxies_output$female_se)
IVW_weights <- F_SHBG_proxies_output$female_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_FEMALE_BETA_CAD ~ allele_matching$ABS_BETA_SHBG- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red", lwd=1.6)


# EGGER
abline(a = 0.004, b = 0.128, col = "blue", lty = 1, lwd=1.6)

# MEDIAN 

legend("topright", legend = c("IWV method", "MR-Egger method"),
       col = c("red", "blue"), lty = c(1, 1), lwd = c(1.6, 1.6))


















############ MALE SHBG ####################################


plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_MALE_BETA_CAD, pch = 16, cex = 0.7,
     xlab = "SNP effect on SHBG",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Male SHBG")

# Add error bars
segments(
  x0 = allele_matching$ABS_BETA_SHBG,
  y0 = allele_matching$HARM_MALE_BETA_CAD - allele_matching$male_se_CAD,
  x1 = allele_matching$ABS_BETA_SHBG,
  y1 = allele_matching$HARM_MALE_BETA_CAD + allele_matching$male_se_CAD,
  col = "black"
)

segments(
  x0 = allele_matching$ABS_BETA_SHBG - allele_matching$SE_SHBG, 
  y0 = allele_matching$HARM_MALE_BETA_CAD,
  x1 = allele_matching$ABS_BETA_SHBG + allele_matching$SE_SHBG, 
  y1 = allele_matching$HARM_MALE_BETA_CAD,
  col = "black"
)

# adding the lines of the different models ############################
# IVW

M_SHBG_proxies_output$male_se <- as.numeric(M_SHBG_proxies_output$male_se)
IVW_weights <- M_SHBG_proxies_output$male_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_SHBG- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red", lwd=1.6)


# EGGER
abline(a = 0.003, b = 0.081, col = "blue", lty = 1, lwd=1.6)

# MEDIAN 

legend("topright", legend = c("IWV method", "MR-Egger method"),
       col = c("red", "blue"), lty = c(1, 1), lwd = c(1.6, 1.6))




snp_x <- 0.12
snp_y <- 0.0074
snp_label <- "rs1799941"

text(snp_x, snp_y, snp_label, col = "blue", pos = 1, cex = 0.7)



###### FEMALE TESTOSTERONE ####################################################


# setting up the basic plot ########################################


plot(allele_matching$ABS_BETA_T, allele_matching$HARM_FEMALE_BETA_CAD, pch = 16, cex = 0.7,
     xlab = "SNP effect on Testosterone",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Female testosterone")

# Add error bars
segments(
  x0 = allele_matching$ABS_BETA_T,
  y0 = allele_matching$HARM_FEMALE_BETA_CAD - allele_matching$female_se_CAD,
  x1 = allele_matching$ABS_BETA_T,
  y1 = allele_matching$HARM_FEMALE_BETA_CAD + allele_matching$female_se_CAD,
  col = "black"
)

segments(
  x0 = allele_matching$ABS_BETA_T - allele_matching$SE_T, 
  y0 = allele_matching$HARM_FEMALE_BETA_CAD,
  x1 = allele_matching$ABS_BETA_T + allele_matching$SE_T,
  y1 = allele_matching$HARM_FEMALE_BETA_CAD,
  col = "black"
)

# adding the lines of the different models ############################
# IVW

F_T_proxies_output$female_se <- as.numeric(F_T_proxies_output$female_se)
IVW_weights <- F_T_proxies_output$female_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_FEMALE_BETA_CAD ~ allele_matching$ABS_BETA_T- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red", lwd=1.6)


# EGGER
abline(a = 0.002, b = -0.058, col = "blue", lty = 1, lwd=1.6)

# MEDIAN 

legend("topright", legend = c("IWV method", "MR-Egger method"),
       col = c("red", "blue"), lty = c(1, 1), lwd = c(1.6, 1.6))












###### MALE TESTOSTERONE ####################################################


# setting up the basic plot ########################################


plot(allele_matching$ABS_BETA_T, allele_matching$HARM_MALE_BETA_CAD, pch = 16, cex = 0.7,
     xlab = "SNP effect on Testosterone",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Male testosterone")

# Add error bars
segments(
  x0 = allele_matching$ABS_BETA_T,
  y0 = allele_matching$HARM_MALE_BETA_CAD - allele_matching$male_se_CAD,
  x1 = allele_matching$ABS_BETA_T,
  y1 = allele_matching$HARM_MALE_BETA_CAD + allele_matching$male_se_CAD,
  col = "black"
)

segments(
  x0 = allele_matching$ABS_BETA_T - allele_matching$SE_T, 
  y0 = allele_matching$HARM_MALE_BETA_CAD,
  x1 = allele_matching$ABS_BETA_T + allele_matching$SE_T,
  y1 = allele_matching$HARM_MALE_BETA_CAD,
  col = "black"
)

# adding the lines of the different models ############################
# IVW

M_T_proxies_output$male_se <- as.numeric(M_T_proxies_output$male_se)
IVW_weights <- M_T_proxies_output$male_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_T- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red", lwd=1.6)


# EGGER
abline(a = 0.002, b = -0.213, col = "blue", lty = 1, lwd=1.6)

# MEDIAN 

legend("topright", legend = c("IWV method", "MR-Egger method"),
       col = c("red", "blue"), lty = c(1, 1), lwd = c(1.6, 1.6))

