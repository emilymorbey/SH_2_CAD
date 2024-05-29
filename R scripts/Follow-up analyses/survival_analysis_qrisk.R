library(survival)
library(lubridate)
library(tidyverse)
surv <- read.csv("data_participant_surv_2.csv")

colnames(surv) <- c("IID", "T", "CAD", "AGERECRUIT", "MONTHBIRTH", "YEARBIRTH", "LTF", "DATEASSESSMENT")
surv$CADBIN <- ifelse(surv$CAD == "" ,0, 1)
surv$LTFBIN <- ifelse(surv$LTF == "" ,0,1)

# adding censoring date
surv$censdate <- Sys.Date()
surv$censdate[surv$LTFBIN == 1] <- surv$LTF[surv$LTFBIN == 1]



# changing the CAD date into a month and a year so we can match it with birth
# date 


# Convert to Date object
surv$CAD <- ymd(surv$CAD)

# Extract month and day
surv$cadmonth <- month(surv$CAD)
surv$cadday <- day(surv$CAD)
surv$cadyear <- year(surv$CAD)

surv$censyear <- year(surv$censdate)


# changing months to numbers for month of birth
surv$monthbirthnum <- as.integer(factor(surv$MONTHBIRTH, levels = month.name))



# creating a censoring variable 
surv$censored <- ifelse(surv$CADBIN == 0, 1, 0)

# creating a year of recruitment variable - this will be inaccurate 
surv$year_of_recruitment <- surv$YEARBIRTH + surv$AGERECRUIT

# creating a time to event variable
surv$timetoCAD <- surv$cadyear-surv$YEARBIRTH

# creating a time to censoring variable 
surv$timetoCENSOR <- ifelse(surv$censored == 1, surv$censyear - surv$YEARBIRTH, NA)

# creating a general time to event variable 
surv$timetoEVENT <- ifelse(is.na(surv$timetoCENSOR), surv$timetoCAD, surv$timetoCENSOR)


# create a variable for testosterone deficiency 
surv$testosterone_deficiency <- ifelse(surv$T < 12, 1, 0)


# need complete cases for testosterone deficiency 
surv <- surv[complete.cases(surv$T), ]






# coxph

CADsurv <- Surv(surv$timetoEVENT, surv$CADBIN)

surv$testosterone_deficiency <- factor(surv$testosterone_deficiency, levels=c(0,1))
km1 <- survfit(CADsurv~surv$testosterone_deficiency)
summary(km1, censored = T)
plot(km1, lty = c(1, 2), col = c("blue", "red"),
     xlab = "Time", ylab = "Survival Probability",
     main = "Kaplan-Meier Curves for Testosterone Deficiency")
legend("bottomright", legend = c("Non-Testosterone Deficient", "Testosterone Deficient"),
       col = c("blue", "red"), lty = c(1, 2))


plot(km1, lty = c(1, 2), col = c("blue", "red"),
     xlab = "Time", ylab = "Survival Probability",
     main = "Kaplan-Meier Curves for Testosterone Deficiency",
     lwd = 2, cex.lab = 1.5)+ # Adjust line thickness and axis label size
legend("bottomright", legend = c("Non-Testosterone Deficient", "Testosterone Deficient"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)

cox1 <- coxph(CADsurv~surv$testosterone_deficiency+surv$AGERECRUIT)
summary(cox1)

# Plot Kaplan-Meier curves using ggsurvplot




















################################################################################
##################### ADDING IN COVARIATES #####################################
################################################################################

surv <- read.csv("data_participant_surv_2.csv")




########### CALCULATING SURVIVOR VARIABLES - DATES ETC. #########################

colnames(surv) <- c("IID", "T", "CAD", "AGERECRUIT", "MONTHBIRTH", "YEARBIRTH", "LTF", "DATEASSESSMENT")
surv$CADBIN <- ifelse(surv$CAD == "" ,0, 1)
surv$LTFBIN <- ifelse(surv$LTF == "" ,0,1)

# adding censoring date
surv$censdate <- Sys.Date()
surv$censdate[surv$LTFBIN == 1] <- surv$LTF[surv$LTFBIN == 1]

# changing the CAD date into a month and a year so we can match it with birth
# date 
# Convert to Date object
surv$CAD <- ymd(surv$CAD)

# Extract month and day
surv$cadmonth <- month(surv$CAD)
surv$cadday <- day(surv$CAD)
surv$cadyear <- year(surv$CAD)
surv$censyear <- year(surv$censdate)

# changing months to numbers for month of birth
surv$monthbirthnum <- as.integer(factor(surv$MONTHBIRTH, levels = month.name))

# creating a censoring variable 
surv$censored <- ifelse(surv$CADBIN == 0, 1, 0)

# creating a year of recruitment variable - this will be inaccurate 
surv$year_of_recruitment <- surv$YEARBIRTH + surv$AGERECRUIT

# creating a time to event variable
surv$timetoCAD <- surv$cadyear-surv$YEARBIRTH

# creating a time to censoring variable 
surv$timetoCENSOR <- ifelse(surv$censored == 1, surv$censyear - surv$YEARBIRTH, NA)

# creating a general time to event variable 
surv$timetoEVENT <- ifelse(is.na(surv$timetoCENSOR), surv$timetoCAD, surv$timetoCENSOR)

# create a variable for testosterone deficiency 
surv$testosterone_deficiency <- ifelse(surv$T < 12, 1, 0)

# need complete cases for testosterone deficiency 
surv <- surv[complete.cases(surv$T), ]






######################### MERGING WITH COVARIATES #############################

qrisk1 <- read.csv("qriskdata/qriskpt1.csv")
qrisk2 <- read.csv("qriskdata/qriskpt2.csv")

names(qrisk1)[names(qrisk1) == "Participant.ID"] <- "IID"
names(qrisk2)[names(qrisk2) == "Participant.ID"] <- "IID"

qrisksurv <- merge(merge(qrisk1, qrisk2, by = "IID", all = TRUE), surv, by = "IID", all = TRUE)

colnames(qrisksurv) <- c("IID", "Agerecruit", "Deprivation", "Ethnicity", 
                         "Cholesterol_aliquot", "Cholesterol_level", 
                         "HDL_instance1", "HDL_instance0", "SBP", 
                         "Diabetes_diagnosed_doctor", "Date_Type1_D", 
                         "Date_Type2_D", "Ever_smoked", "Current_smoking", 
                         "Num_cigarettes_daily", "Cholesterol_medication_0", 
                         "Cholesterol_medication_1", "Microalbumin_inst_0", 
                         "Microalbumin_inst_1", "Creatinine_inst_0", "Ill_fath", 
                         "Fath_alive", "Ill_moth", "Ill_sib", "non-canc_illness", 
                         "Date_a_fib", "BMI1", "BMI2", "DATE_SLE", "DATE_MIGRAINE", 
                         "DATE_SEROPOS_ARTHRITIS", "DATE_OTHER_ARTHRITIS", 
                         "DATE_JUV_ARTHRITIS", "DATE_SCHIZO", "DATE_BIPOLAR", 
                         "DATE_DEPRESSIVE_RECUR", "DATE_DEPRESSIVE_EP", 
                         "TREATMENT_MEDS", "NON_CANC_ILLNESS", "T", "CAD", 
                         "AGERECRUIT", "MONTHBIRTH", "YEARBIRTH", "LTF", 
                         "DATEASSESSMENT", "CADBIN", "LTFBIN", "CENSDATE", 
                         "CADMONTH", "CADDAY", "CADYEAR", "CENSYEAR", 
                         "MONTHBIRTHNUM", "CENSORED", "YEAR_RECRUIT", 
                         "timetoCAD", "timetoCENSOR", "timetoEVENT", 
                         "TESTOSTERONE_DEFICIENCY")



# removing people with no info on age or entrance to study 
qrisksurv <- qrisksurv[!is.na(qrisksurv$AGERECRUIT), ]


# MAKING NEW VARIABLES TO SUMMARISE THE DIFFERENT COVARIATES 

qrisksurv <- qrisksurv[, -which(names(qrisksurv) %in% c("Diabetes_diagnosed_doctor", "Microalbumin_inst_1", "Fath_alive"))]

# diabetes 


qrisksurv$TYPE1DBIN <- ifelse(qrisksurv$Date_Type1_D == "" ,0,1)
qrisksurv$TYPE2DBIN <- ifelse(qrisksurv$Date_Type2_D == "" ,0,1)
qrisksurv$TYPE1DBIN <- factor(qrisksurv$TYPE1DBIN, levels=c(0,1))
qrisksurv$TYPE2DBIN <- factor(qrisksurv$TYPE2DBIN, levels=c(0,1))

qrisksurv <- qrisksurv[, -which(names(qrisksurv) %in% c("Date_Type1_D", "Date_Type2_D", "Cholesterol_aliquot", "HDL_instance1"))]

table(qrisksurv$TYPE1DBIN)

# cholesterol medication 
table(qrisksurv$Cholesterol_medication_0)

qrisksurv$CholesterolBIN <- ifelse(qrisksurv$Cholesterol_medication_0 == "Cholesterol loweing medication" | 
                                      qrisksurv$Cholesterol_medication_0 =="Cholesterol lowering medication|Blood pressure medication|Insulin" |
                                      qrisksurv$Cholesterol_medication_0 == "Cholesterol lowering medication|Blood pressure medication" | 
                                      qrisksurv$Cholesterol_medication_0 == "Cholesterol lowering medication|Insulin" ,1,0)

qrisksurv$CholesterolBIN <- factor(qrisksurv$CholesterolBIN, levels=c(0,1))


# atrial fibrilation
table(qrisksurv$`non-canc_illness`)

qrisksurv$AFIBBIN1 <- ifelse(qrisksurv$`non-canc_illness` == "atrial fibrillation" ,1,0)
qrisksurv$AFIBBIN1 <- factor(qrisksurv$AFIBBIN1, levels=c(0,1))

# systemic lupus erythematosis 
qrisksurv$SLEBIN <- ifelse(qrisksurv$`non-canc_illness` == "systemic lupus erythematosis/sle" ,1,0)
qrisksurv$SLEBIN <- factor(qrisksurv$SLEBIN, levels=c(0,1))

# migraine
qrisksurv$MIGRAINEBIN1 <- ifelse(qrisksurv$`non-canc_illness` == "migraine" ,1,0)
qrisksurv$MIGRAINEBIN1 <- factor(qrisksurv$MIGRAINEBIN1, levels=c(0,1))

# mental illnesses 
qrisksurv$MENTALBIN <- ifelse(qrisksurv$`non-canc_illness` %in% c("depression", "deliberate self-harm/suicide attempt", "schizophrenia"), 1, 0)
qrisksurv$MENTALBIN <- factor(qrisksurv$MENTALBIN, levels=c(0,1))

# erectile dysfunction 
qrisksurv$EDBIN <- ifelse(qrisksurv$`non-canc_illness` == "erectile dysfunction / impotence" ,1,0)
qrisksurv$EDBIN <- factor(qrisksurv$EDBIN, levels=c(0,1))

# atypical antipsychotic medication 
qrisksurv$ANTIPSYCHOTICMEDBIN <- ifelse(qrisksurv$TREATMENT_MEDS %in% c("risperidone", 
                                        "quetiapine", "olanzapine", "aripiprazole", 
                                        "clozapine", "haloperidol") ,1,0)
qrisksurv$ANTIPSYCHOTICMEDBIN <- factor(qrisksurv$ANTIPSYCHOTICMEDBIN, levels=c(0,1))

# steroids  
qrisksurv$STEROIDSBIN <- ifelse(qrisksurv$TREATMENT_MEDS %in% c("prednisolone", "betamethasone",
                                                                "dexamethasone", "hydrocortisone") ,1,0)
qrisksurv$STEROIDSBIN <- factor(qrisksurv$STEROIDSBIN, levels=c(0,1))

# smoking - this one is more tricky 

table(qrisksurv$Ever_smoked)

qrisksurv$exsmoker <- ifelse(qrisksurv$Ever_smoked == "Yes" & qrisksurv$Current_smoking=="No" ,"ex-smoker",NA)
qrisksurv$nonsmoker <- ifelse(qrisksurv$Ever_smoked == "No" ,"non-smoker",NA)


qrisksurv$SmokingCategory <- ifelse(qrisksurv$Num_cigarettes_daily < 10 & qrisksurv$Num_cigarettes_daily > 0, "Light Smoker",
                                    ifelse(qrisksurv$Num_cigarettes_daily >= 10 & qrisksurv$Num_cigarettes_daily < 20, "Moderate Smoker", 
                                           ifelse(qrisksurv$Num_cigarettes_daily >= 20, "Heavy Smoker", NA)))

table(qrisksurv$SmokingCategory)


qrisksurv$UKBBSMOKING <- qrisksurv %>%
  mutate_at(vars(starts_with("Col")), ~na_if(., "0")) %>%
  mutate(Consolidated = coalesce(nonsmoker, exsmoker, SmokingCategory)) %>%
  select(Consolidated)

qrisksurv$UKBBSMOKING <- factor(qrisksurv$UKBBSMOKING, levels=c("non-smoker","ex-smoker", "Light Smoker", "Moderate Smoker", "Heavy Smoker"))


qrisksurv[is.na(qrisksurv)] <- "Unknown"


### TRYING TO ADJUST FOR A FEW THINGS

CADsurv <- Surv(qrisksurv$timetoEVENT, qrisksurv$CADBIN)


qrisksurv$TESTOSTERONE_DEFICIENCY <- factor(qrisksurv$TESTOSTERONE_DEFICIENCY, levels=c(0,1))
km1 <- survfit(CADsurv~qrisksurv$TESTOSTERONE_DEFICIENCY + qrisksurv$AGERECRUIT + qrisksurv$Deprivation + qrisksurv$Cholesterol_level +
                 qrisksurv$SBP + qrisksurv$BMI1)


qrisksurv$TYPE1DBIN + qrisksurv$TYPE2DBIN + 
  qrisksurv$CholesterolBIN1 + qrisksurv$AFIBBIN1 + qrisksurv$SLEBIN + qrisksurv$MIGRAINEBIN1 +
  qrisksurv$MENTALBIN + qrisksurv$EDBIN + qrisksurv$ANTIPSYCHOTICMEDBIN + qrisksurv$STEROIDSBIN +
  qrisksurv$UKBBSMOKING
summary(km1, censored = T)
plot(km1, lty = c(1, 2), col = c("blue", "red"),
     xlab = "Time", ylab = "Survival Probability",
     main = "Kaplan-Meier Curves for Testosterone Deficiency")
legend("bottomright", legend = c("Non-Testosterone Deficient", "Testosterone Deficient"),
       col = c("blue", "red"), lty = c(1, 2))









































library(survival)
library(lubridate)
library(tidyverse)
library(dplyr)
surv <- read.csv("/mnt/project/data/data_participant_surv_2.csv")

colnames(surv) <- c("IID", "T", "CAD", "AGERECRUIT", "MONTHBIRTH", "YEARBIRTH", "LTF", "DATEASSESSMENT")
surv$CADBIN <- ifelse(surv$CAD == "" ,0, 1)
surv$LTFBIN <- ifelse(surv$LTF == "" ,0,1)

# adding censoring date
surv$censdate <- Sys.Date()
surv$censdate[surv$LTFBIN == 1] <- surv$LTF[surv$LTFBIN == 1]



# changing the CAD date into a month and a year so we can match it with birth
# date 


# Convert to Date object
surv$CAD <- ymd(surv$CAD)

# Extract month and day
surv$cadmonth <- month(surv$CAD)
surv$cadday <- day(surv$CAD)
surv$cadyear <- year(surv$CAD)

surv$censyear <- year(surv$censdate)


# changing months to numbers for month of birth
surv$monthbirthnum <- as.integer(factor(surv$MONTHBIRTH, levels = month.name))



# creating a censoring variable 
surv$censored <- ifelse(surv$CADBIN == 0, 1, 0)

# creating a year of recruitment variable - this will be inaccurate 
surv$year_of_recruitment <- surv$YEARBIRTH + surv$AGERECRUIT

# creating a time to event variable
surv$timetoCAD <- surv$cadyear-surv$YEARBIRTH

# creating a time to censoring variable 
surv$timetoCENSOR <- ifelse(surv$censored == 1, surv$censyear - surv$YEARBIRTH, NA)

# creating a general time to event variable 
surv$timetoEVENT <- ifelse(is.na(surv$timetoCENSOR), surv$timetoCAD, surv$timetoCENSOR)


# create a variable for testosterone deficiency 
surv$testosterone_deficiency <- ifelse(surv$T < 12, 1, 0)


# need complete cases for testosterone deficiency 
surv <- surv[complete.cases(surv$T), ]






# coxph

CADsurv <- Surv(surv$timetoEVENT, surv$CADBIN)

surv$testosterone_deficiency <- factor(surv$testosterone_deficiency, levels=c(0,1))
km1 <- survfit(CADsurv~surv$testosterone_deficiency)
summary(km1, censored = T)
plot(km1, lty = c(1, 2), col = c("blue", "red"),
     xlab = "Time", ylab = "Survival Probability",
     main = "Kaplan-Meier Curves for Testosterone Deficiency")
legend("bottomright", legend = c("Non-Testosterone Deficient", "Testosterone Deficient"),
       col = c("blue", "red"), lty = c(1, 2))

cox1 <- coxph(CADsurv~surv$testosterone_deficiency+surv$AGERECRUIT)
summary(cox1)

# Plot Kaplan-Meier curves using ggsurvplot




















################################################################################
##################### ADDING IN COVARIATES #####################################
################################################################################

surv <- read.csv("/mnt/project/data/data_participant_surv_2.csv")




########### CALCULATING SURVIVOR VARIABLES - DATES ETC. #########################

colnames(surv) <- c("IID", "T", "CAD", "AGERECRUIT", "MONTHBIRTH", "YEARBIRTH", "LTF", "DATEASSESSMENT")
surv$CADBIN <- ifelse(surv$CAD == "" ,0, 1)
surv$LTFBIN <- ifelse(surv$LTF == "" ,0,1)

# adding censoring date
surv$censdate <- Sys.Date()
surv$censdate[surv$LTFBIN == 1] <- surv$LTF[surv$LTFBIN == 1]

# changing the CAD date into a month and a year so we can match it with birth
# date 
# Convert to Date object
surv$CAD <- ymd(surv$CAD)

# Extract month and day
surv$cadmonth <- month(surv$CAD)
surv$cadday <- day(surv$CAD)
surv$cadyear <- year(surv$CAD)
surv$censyear <- year(surv$censdate)

# changing months to numbers for month of birth
surv$monthbirthnum <- as.integer(factor(surv$MONTHBIRTH, levels = month.name))

# creating a censoring variable 
surv$censored <- ifelse(surv$CADBIN == 0, 1, 0)

# creating a year of recruitment variable - this will be inaccurate 
surv$year_of_recruitment <- surv$YEARBIRTH + surv$AGERECRUIT

# creating a time to event variable
surv$timetoCAD <- surv$cadyear-surv$YEARBIRTH

# creating a time to censoring variable 
surv$timetoCENSOR <- ifelse(surv$censored == 1, surv$censyear - surv$YEARBIRTH, NA)

# creating a general time to event variable 
surv$timetoEVENT <- ifelse(is.na(surv$timetoCENSOR), surv$timetoCAD, surv$timetoCENSOR)

# create a variable for testosterone deficiency 
surv$testosterone_deficiency <- ifelse(surv$T < 12, 1, 0)

# need complete cases for testosterone deficiency 
surv <- surv[complete.cases(surv$T), ]






######################### MERGING WITH COVARIATES #############################

qrisk1 <- read.csv("/mnt/project/data/qriskpt1.csv")
qrisk2 <- read.csv("/mnt/project/data/qriskpt2.csv")

names(qrisk1)[names(qrisk1) == "Participant.ID"] <- "IID"
names(qrisk2)[names(qrisk2) == "Participant.ID"] <- "IID"

qrisksurv <- merge(merge(qrisk1, qrisk2, by = "IID", all = TRUE), surv, by = "IID", all = TRUE)

colnames(qrisksurv) <- c("IID", "Agerecruit", "Deprivation", "Ethnicity", 
                         "Cholesterol_aliquot", "Cholesterol_level", 
                         "HDL_instance1", "HDL_instance0", "SBP", 
                         "Diabetes_diagnosed_doctor", "Date_Type1_D", 
                         "Date_Type2_D", "Ever_smoked", "Current_smoking", 
                         "Num_cigarettes_daily", "Cholesterol_medication", 
                         "Cholesterol_medication_2", "Microalbumin_inst_0", 
                         "Microalbumin_inst_1", "Creatinine_inst_0", "Ill_fath", 
                         "Fath_alive", "Ill_moth", "Ill_sib", "non-canc_illness", 
                         "Date_a_fib", "BMI1", "BMI2", "DATE_SLE", "DATE_MIGRAINE", 
                         "DATE_SEROPOS_ARTHRITIS", "DATE_OTHER_ARTHRITIS", 
                         "DATE_JUV_ARTHRITIS", "DATE_SCHIZO", "DATE_BIPOLAR", 
                         "DATE_DEPRESSIVE_RECUR", "DATE_DEPRESSIVE_EP", 
                         "TREATMENT_MEDS", "NON_CANC_ILLNESS", "T", "CAD", 
                         "AGERECRUIT", "MONTHBIRTH", "YEARBIRTH", "LTF", 
                         "DATEASSESSMENT", "CADBIN", "LTFBIN", "CENSDATE", 
                         "CADMONTH", "CADDAY", "CADYEAR", "CENSYEAR", 
                         "MONTHBIRTHNUM", "CENSORED", "YEAR_RECRUIT", 
                         "timetoCAD", "timetoCENSOR", "timetoEVENT", 
                         "TESTOSTERONE_DEFICIENCY")



# removing people with no info on age or entrance to study 
qrisksurv <- qrisksurv[!is.na(qrisksurv$AGERECRUIT), ]


# MAKING NEW VARIABLES TO SUMMARISE THE DIFFERENT COVARIATES 

qrisksurv <- qrisksurv[, -which(names(qrisksurv) %in% c("Diabetes_diagnosed_doctor", "Microalbumin_inst_1", "Fath_alive"))]

# diabetes 
qrisksurv$TYPE1DBIN <- ifelse(qrisksurv$Date_Type1_D == "" ,0,1)
qrisksurv$TYPE2DBIN <- ifelse(qrisksurv$Date_Type2_D == "" ,0,1)
qrisksurv$TYPE1DBIN <- factor(qrisksurv$TYPE1DBIN, levels=c(0,1))
qrisksurv$TYPE2DBIN <- factor(qrisksurv$TYPE2DBIN, levels=c(0,1))

qrisksurv <- qrisksurv[, -which(names(qrisksurv) %in% c("Date_Type1_D", "Date_Type2_D", "Cholesterol_aliquot", "HDL_instance1"))]


# cholesterol medication 
qrisksurv$CholesterolBIN1 <- ifelse(qrisksurv$Cholesterol_medication == "None of the above" ,0,1)
qrisksurv$CholesterolBIN1 <- factor(qrisksurv$CholesterolBIN1, levels=c(0,1))

# atrial fibrilation
qrisksurv$AFIBBIN1 <- ifelse(qrisksurv$`non-canc_illness` == "atrial fibrillation" ,1,0)
qrisksurv$AFIBBIN1 <- factor(qrisksurv$AFIBBIN1, levels=c(0,1))

# systemic lupus erythematosis 
qrisksurv$SLEBIN <- ifelse(qrisksurv$`non-canc_illness` == "systemic lupus erythematosis/sle" ,1,0)
qrisksurv$SLEBIN <- factor(qrisksurv$SLEBIN, levels=c(0,1))

# migraine
qrisksurv$MIGRAINEBIN1 <- ifelse(qrisksurv$`non-canc_illness` == "migraine" ,1,0)
qrisksurv$MIGRAINEBIN1 <- factor(qrisksurv$MIGRAINEBIN1, levels=c(0,1))

# mental illnesses 
qrisksurv$MENTALBIN <- ifelse(qrisksurv$`non-canc_illness` %in% c("depression", "deliberate self-harm/suicide attempt", "schizophrenia"), 1, 0)
qrisksurv$MENTALBIN <- factor(qrisksurv$MENTALBIN, levels=c(0,1))

# erectile dysfunction 
qrisksurv$EDBIN <- ifelse(qrisksurv$`non-canc_illness` == "erectile dysfunction / impotence" ,1,0)
qrisksurv$EDBIN <- factor(qrisksurv$EDBIN, levels=c(0,1))

# atypical antipsychotic medication 
qrisksurv$ANTIPSYCHOTICMEDBIN <- ifelse(qrisksurv$TREATMENT_MEDS %in% c("risperidone", 
                                                                        "quetiapine", "olanzapine", "aripiprazole", 
                                                                        "clozapine", "haloperidol") ,1,0)
qrisksurv$ANTIPSYCHOTICMEDBIN <- factor(qrisksurv$ANTIPSYCHOTICMEDBIN, levels=c(0,1))

# steroids  
qrisksurv$STEROIDSBIN <- ifelse(qrisksurv$TREATMENT_MEDS %in% c("prednisolone", "betamethasone",
                                                                "dexamethasone", "hydrocortisone") ,1,0)
qrisksurv$STEROIDSBIN <- factor(qrisksurv$STEROIDSBIN, levels=c(0,1))


# arthritis 
qrisksurv$ARTHBIN <- ifelse(qrisksurv$DATE_SEROPOS_ARTHRITIS == "" ,0,1)
qrisksurv$ARTHBIN <- factor(qrisksurv$ARTHBIN, levels=c(0,1))


# bp treatment

sum(qrisksurv$TREATMENT_MEDS == "olmesartan")

qrisksurv$BPTREATBIN <- ifelse(qrisksurv$TREATMENT_MEDS %in% c("enalapril", "lisinopril",
                                                               "perindopril", "ramipril", 
                                                               "losartan", "valsartan", "olmesartan", 
                                                               "amlodipine", "felodipine", "nifedipine", "diltiazem", 
                                                               "verapami", "indapamide", "bendroflumethiazide", "atenolol", 
                                                               "bisoprolol") ,1,0)

qrisksurv$BPTREATBIN <- factor(qrisksurv$STEROIDSBIN, levels=c(0,1))



# ethnicity 




# smoking - this one is more tricky 

qrisksurv$exsmoker <- ifelse(qrisksurv$Ever_smoked == "Yes" & qrisksurv$Current_smoking=="No" ,"2",NA)
qrisksurv$nonsmoker <- ifelse(qrisksurv$Ever_smoked == "No" ,"1",NA)

qrisksurv$Num_cigarettes_daily <- as.numeric(qrisksurv$Num_cigarettes_daily)

qrisksurv$SmokingCategory <- ifelse(qrisksurv$Num_cigarettes_daily < 10 & qrisksurv$Num_cigarettes_daily > 1, "3",
                                    ifelse(qrisksurv$Num_cigarettes_daily >= 10 & qrisksurv$Num_cigarettes_daily < 20, "4", 
                                           ifelse(qrisksurv$Num_cigarettes_daily >= 20, "5", NA)))



# Remove NA values and replace with empty strings
qrisksurv$exsmoker[is.na(qrisksurv$exsmoker)] <- ""
qrisksurv$nonsmoker[is.na(qrisksurv$nonsmoker)] <- ""
qrisksurv$SmokingCategory[is.na(qrisksurv$SmokingCategory)] <- ""

# Combine columns into UKBBSMOKING with no white space
qrisksurv$UKBBSMOKING <- paste0(qrisksurv$exsmoker, qrisksurv$nonsmoker, qrisksurv$SmokingCategory)


qrisksurv$UKBBSMOKING <- as.factor(qrisksurv$UKBBSMOKING)

qrisksurv %>%
  mutate(UKBBSMOKING = factor(UKBBSMOKING,
                              levels = c("","1","2","3","4","5")))


# gathering the relevant columns



qrisksurv <- select(qrisksurv, 
                    IID, AGERECRUIT, CADBIN, timetoEVENT, TESTOSTERONE_DEFICIENCY, 
                    BMI1, TYPE1DBIN, TYPE2DBIN, EDBIN, SLEBIN, AFIBBIN1, UKBBSMOKING, 
                    MENTALBIN, MIGRAINEBIN1, STEROIDSBIN, ANTIPSYCHOTICMEDBIN, 
                    CholesterolBIN1)




### ADJUSTING FOR AGE AND ERECTILE DYSFUNCTION 

CADsurv <- Surv(qrisksurv$timetoEVENT, qrisksurv$CADBIN)


qrisksurv$TESTOSTERONE_DEFICIENCY <- factor(qrisksurv$TESTOSTERONE_DEFICIENCY, levels=c(0,1))
km1 <- survfit(CADsurv~qrisksurv$TESTOSTERONE_DEFICIENCY + qrisksurv$AGERECRUIT + qrisksurv$EDBIN)


summary(km1, censored = T)

plot(km1, lty = c(1, 2), col = c("blue", "red"),
     xlab = "Time", ylab = "Survival Probability",
     main = "Kaplan-Meier Curves for Testosterone Deficiency")
legend("bottomright", legend = c("Non-Testosterone Deficient", "Testosterone Deficient"),
       col = c("blue", "red"), lty = c(1, 2))




coxED <- coxph(CADsurv~qrisksurv$TESTOSTERONE_DEFICIENCY+qrisksurv$AGERECRUIT+qrisksurv$EDBIN)
summary(coxED)







# ADJUSTING FOR AGE, ED AND DIABETES

CADsurv <- Surv(qrisksurv$timetoEVENT, qrisksurv$CADBIN)


qrisksurv$TESTOSTERONE_DEFICIENCY <- factor(qrisksurv$TESTOSTERONE_DEFICIENCY, levels=c(0,1))
km1 <- survfit(CADsurv~qrisksurv$TESTOSTERONE_DEFICIENCY + qrisksurv$AGERECRUIT + qrisksurv$EDBIN +
                 qrisksurv$TYPE1DBIN + qrisksurv$TYPE2DBIN)


coxEDD <- coxph(CADsurv~qrisksurv$TESTOSTERONE_DEFICIENCY+qrisksurv$AGERECRUIT+qrisksurv$EDBIN+
                  qrisksurv$TYPE1DBIN + qrisksurv$TYPE2DBIN)
summary(coxEDD)







# ADJUSTING FOR AGE, ED, DIABETES, CHOLESTEROL AND AFIB

CADsurv <- Surv(qrisksurv$timetoEVENT, qrisksurv$CADBIN)


qrisksurv$TESTOSTERONE_DEFICIENCY <- factor(qrisksurv$TESTOSTERONE_DEFICIENCY, levels=c(0,1))
km1 <- survfit(CADsurv~qrisksurv$TESTOSTERONE_DEFICIENCY + qrisksurv$AGERECRUIT + qrisksurv$EDBIN +
                 qrisksurv$TYPE1DBIN + qrisksurv$TYPE2DBIN + qrisksurv$CholesterolBIN1 + qrisksurv$AFIBBIN1)


coxEDDCA <- coxph(CADsurv~qrisksurv$TESTOSTERONE_DEFICIENCY+qrisksurv$AGERECRUIT+qrisksurv$EDBIN+
                    qrisksurv$TYPE1DBIN + qrisksurv$TYPE2DBIN + qrisksurv$CholesterolBIN1 + qrisksurv$AFIBBIN1)
summary(coxEDDCA)





# ADJUSTING FOR AGE, ED, DIABETES, CHOLESTEROL, AFIB 

CADsurv <- Surv(qrisksurv$timetoEVENT, qrisksurv$CADBIN)


qrisksurv$TESTOSTERONE_DEFICIENCY <- factor(qrisksurv$TESTOSTERONE_DEFICIENCY, levels=c(0,1))
km1 <- survfit(CADsurv~qrisksurv$TESTOSTERONE_DEFICIENCY + qrisksurv$AGERECRUIT + qrisksurv$EDBIN +
                 qrisksurv$TYPE1DBIN + qrisksurv$TYPE2DBIN + qrisksurv$CholesterolBIN1 + qrisksurv$AFIBBIN1
               + qrisksurv$SLEBIN + qrisksurv$MIGRAINEBIN1)


coxEDDCASM <- coxph(CADsurv~qrisksurv$TESTOSTERONE_DEFICIENCY+qrisksurv$AGERECRUIT+qrisksurv$EDBIN+
                      qrisksurv$TYPE1DBIN + qrisksurv$TYPE2DBIN + qrisksurv$CholesterolBIN1 + qrisksurv$AFIBBIN1
                    + qrisksurv$SLEBIN + qrisksurv$MIGRAINEBIN1)
summary(coxEDDCASM)




# ADJUSTING FOR ALL 

CADsurv <- Surv(qrisksurv$timetoEVENT, qrisksurv$CADBIN)


qrisksurv$TESTOSTERONE_DEFICIENCY <- factor(qrisksurv$TESTOSTERONE_DEFICIENCY, levels=c(0,1))
km1 <- survfit(CADsurv~qrisksurv$TESTOSTERONE_DEFICIENCY + qrisksurv$AGERECRUIT + qrisksurv$EDBIN +
                 qrisksurv$TYPE1DBIN + qrisksurv$TYPE2DBIN + qrisksurv$CholesterolBIN1 + qrisksurv$AFIBBIN1
               + qrisksurv$SLEBIN + qrisksurv$MIGRAINEBIN1 + qrisksurv$MENTALBIN + qrisksurv$EDBIN + qrisksurv$STEROIDSBIN +
                 qrisksurv$UKBBSMOKING)


coxALL <- coxph(CADsurv~qrisksurv$TESTOSTERONE_DEFICIENCY+qrisksurv$AGERECRUIT+qrisksurv$EDBIN+
                  qrisksurv$TYPE1DBIN + qrisksurv$TYPE2DBIN + qrisksurv$CholesterolBIN1 + qrisksurv$AFIBBIN1
                + qrisksurv$SLEBIN + qrisksurv$MIGRAINEBIN1 +qrisksurv$MENTALBIN +
                  qrisksurv$UKBBSMOKING + qrisksurv$STEROIDSBIN + qrisksurv$ANTIPSYCHOTICMEDBIN)



summary(coxALL)


