

library(tidyverse)
library(dplyr)

PRS <- read.table("PRS_MALE_TEST.txt.profile", header = TRUE)


# Load necessary libraries
library(ggplot2)

# Plot histogram of PRS
ggplot(PRS, aes(x = SCORE)) +
  geom_histogram(binwidth = 0.00001, fill = "skyblue", color = "blue") +
  labs(title = "Distribution of Polygenic Risk Scores",
       x = "Polygenic Risk Score",
       y = "Frequency") +
  theme_minimal()






# TRYING TO DO THE ONE SAMPLE MR 

scores <- read.table("PRS_MALE_TEST.txt.profile", header = TRUE)
phenotypes <- read.csv("new_cad_1.csv")
phenotypes2 <- read.csv("new_cad_2.csv")

phenotypes_collapsed <- phenotypes %>%
  group_by(eid) %>%
  summarize(
    diag_icd10 = paste(diag_icd10, collapse = ", "),
    diag_icd9 = paste(diag_icd9, collapse = "")
  ) %>%
  ungroup()

phenotypes_collapsed <- phenotypes_collapsed %>%
  mutate(CAD_ICD10 = if_else(grepl("I21|I210|I211|I212|I213|
                            I214|I219|I21X|I22|I220|I221|I228|I229|
                            I23|I23.1|I23.2|I23.3|I23.4|I23.5|
                            I23.6|I238|I24|I240|I241|I248|I249|
                            I252", phenotypes_collapsed$diag_icd10), 1, 0))


phenotypes_collapsed <- phenotypes_collapsed %>%
  mutate("CAD_ICD9" = if_else(grepl("410|4109|411|4119|
                                    412|4129", phenotypes_collapsed$diag_icd9), 1, 0))


phenotypes_all <- merge(phenotypes_collapsed, phenotypes2, by = "eid")



phenotypes_full <- merge(phenotypes, phenotypes2, by = "eid", all = TRUE)



names(phenotypes_all)[names(phenotypes_all) == "X30850.0.0"] <- "T"
names(phenotypes_all)[names(phenotypes_all) == "X41272.0.0"] <- "OPS"
names(phenotypes_all)[names(phenotypes_all) == "X20002.0.0"] <- "Self-report"
names(phenotypes_all)[names(phenotypes_all) == "eid"] <- "IID"


phenotypes_all <- phenotypes_all %>%
  mutate("CAD_OP" = if_else(grepl("K40|K401|K402|K403|K404|K408|K409|K41|K411|
  K412|K413|K414|K418|K419|K42|K421|K422|K423|K424|K428|K429|K43|K431|K432|K433|
  K434|K438|K439|K44|K441|K442|K448|K449|K45|K451|K452|K453|K454|K455|K456|K458|
  K459|K46|K461|K462|K463|K464|K465|K468|K469|K49|K491|K492|K49.3|K494|K498|K499|
  K50|K501|K502|K50.3|K504|K508|K509|K75|K75.1|K75.2|K75.3|K75.4|K75.8|K75.9", phenotypes_all$OPS), 1, 0))


sum(phenotypes_all$CAD_OP==1)


phenotypes_all <- phenotypes_all %>%
  mutate("CAD_SELFREP" = if_else(grepl("1075|1070|1095|1523", phenotypes_all$`Self-report`), 1, 0))

sum(phenotypes_all$CAD_SELFREP==1)

phenotypes_all$CAD_ICD10 <- as.numeric(phenotypes_all$CAD_ICD10)
phenotypes_all$CAD_ICD9 <- as.numeric(phenotypes_all$CAD_ICD9)
phenotypes_all$CAD_OP <- as.numeric(phenotypes_all$CAD_OP)
phenotypes_all$CAD_SELFREP <- as.numeric(phenotypes_all$CAD_SELFREP)

phenotypes_all$CADBIN <- as.numeric(rowSums(phenotypes_all[, c("CAD_ICD10", "CAD_ICD9", "CAD_OP", "CAD_SELFREP")]) > 0)

sum(phenotypes_all$CADBIN==1)

phenoscores <- merge(phenotypes_all, scores, by = "IID", all = TRUE)


sum(phenoscores$CADBIN=="1", na.rm = TRUE)

# one sample MR
CADBIN <- factor(phenoscores$CADBIN, levels = c("0", "1"), ordered = TRUE)

plot(phenoscores$SCORE, phenoscores$T)

library(ggplot2)
ggplot(phenoscores, aes(x = SCORE, y = T)) +
  geom_point() +  # Add scatterplot points
  geom_smooth(method = "lm", se = FALSE) +  # Add linear model line without confidence interval
  labs(x = "SCORE", y = "T")


# regressing testosterone on the score 
lm_model <- lm(T ~ SCORE, data = phenoscores)
summary(lm_model)

# regressing CAD on the score (logistic)
lm_model_cad <- glm(CADBIN ~ SCORE, data = phenoscores, family = binomial)
summary(lm_model_cad)

# Box plot
ggplot(phenoscores, aes(x = factor(CADBIN), y = SCORE)) +
  geom_boxplot() +
  labs(x = "CADBIN", y = "Genetic score", title = "Box Plot of CAD and polygenic score")

phenoscores$SCORE <- as.numeric(phenoscores$SCORE)

ggplot(phenoscores, aes(x = SCORE, fill = CADBIN)) +
  geom_histogram(binwidth = 1, position = "dodge") +
  facet_wrap(~CADBIN, scales = "free") +
  labs(title = "Histogram of SCORE by CADBIN",
       x = "SCORE",
       y = "Frequency")

# regressing cad on t (logistic)
lm_model_cad_T <- glm(CADBIN ~ T, data = phenoscores, family = binomial)
summary(lm_model_cad_T)

# Box plot
ggplot(phenoscores, aes(x = factor(CADBIN), y = T)) +
  geom_boxplot() +
  labs(x = "CADBIN", y = "Testosterone Levels", title = "Box Plot of CADBIN and Testosterone Levels")





library(ivreg)

iv_mod <- ivreg(CADBIN ~ T | SCORE, data=phenoscores)
summary(iv_mod)








# ONE SAMPLE MR OF BLOOD PRESSURE 

scores <- read.table("PRS_MALE_TEST.txt.profile", header = TRUE)
phenotypes <- read.csv("data_participant_bp.csv")

colnames(phenotypes) <- c("IID", "T", "CAD", "AGE", "BP")


phenoscores <- merge(phenotypes, scores, by = "IID", all = TRUE)
phenoscores <- na.omit(phenoscores)



# one sample MR

plot(phenoscores$SCORE, phenoscores$T)

library(ggplot2)
ggplot(phenoscores, aes(x = SCORE, y = T)) +
  geom_point() +  # Add scatterplot points
  geom_smooth(method = "lm", se = FALSE) +  # Add linear model line without confidence interval
  labs(x = "SCORE", y = "T")

lm_model <- lm(T ~ SCORE, data = phenoscores)
lm_model_BP <- lm(BP ~ SCORE, data = phenoscores)
summary(lm_model_BP)

lm_model_BP_T <- lm(BP ~ T, data = phenoscores)
summary(lm_model_BP_T)

library(ivreg)

iv_mod <- ivreg(BP ~ T | SCORE, data=phenoscores)
summary(iv_mod)




# ADDING AGE INTO THE BP MODEL 
scores <- read.table("PRS_MALE_TEST.txt.profile", header = TRUE)
phenotypes <- read.csv("data_participant_bp.csv")

colnames(phenotypes) <- c("IID", "T", "CAD", "AGE", "BP")

phenotypes$CAD <- as.character(phenotypes$CAD)

phenotypes$CADBIN <- ifelse(phenotypes$CAD == "" ,0, 1)

phenoscores <- merge(phenotypes, scores, by = "IID", all = TRUE)
phenoscores <- na.omit(phenoscores)




# regressing testosterone on the score 
lm_model <- lm(T ~ SCORE + AGE, data = phenoscores)
summary(lm_model)

# regressing CAD on the score (logistic)
lm_model_cad <- glm(CADBIN ~ SCORE + AGE, data = phenoscores, family = binomial)
summary(lm_model_cad)

# regressing cad on t (logistic)
lm_model_cad_T <- glm(CADBIN ~ T + AGE, data = phenoscores, family = binomial)
summary(lm_model_cad_T)


library(ivreg)

iv_mod <- ivreg(CADBIN ~ T+AGE | SCORE+AGE, data=phenoscores)
summary(iv_mod)




# Box plot
ggplot(phenoscores, aes(x = factor(CADBIN), y = T)) +
  geom_boxplot() +
  labs(x = "CADBIN", y = "Testosterone Levels", title = "Box Plot of CADBIN and Testosterone Levels")












