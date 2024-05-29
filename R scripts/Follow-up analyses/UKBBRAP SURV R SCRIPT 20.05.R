
install.packages("data.table")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("ggforce")
install.packages("tableone")
install.packages("survminer")
install.packages("MatchIt")
install.packages("devtools")
install.packages("cli")


library(tidyverse)
library(dplyr)
library(survival)
library(lubridate)
library(stringr)
library(cli)
library(data.table)



################################################################################
##################### ADDING IN COVARIATES #####################################
################################################################################


surv <- read.csv("data_participant_surv_2.csv")

colnames(surv) <- c("IID", "T", "CAD", "AGERECRUIT", "MONTHBIRTH", "YEARBIRTH", "LTF", "DATEASSESSMENT")



################################################################################
######################### PHENOTYPING CAD ######################################
################################################################################


icd10_1 <- read.csv("male_icd10_1.csv")
icd10_2 <- read.csv("male_icd10_2.csv")


patterns <- c("I21", "I210", "I211", "I212", "I213", "I214", "I219", "I21X",
              "I22", "I220", "I221", "I228", "I229", "I23", "I23.1", "I23.2",
              "I23.3", "I23.4", "I23.5", "I23.6", "I238", "I24", "I240",
              "I241", "I248", "I249", "I252")


relevant_icd10 <- icd10_2 %>%
  filter(grepl(paste(patterns, collapse = "|"), diag_icd10))


icd10_dates <- merge(icd10_1, relevant_icd10, by = "eid")
icd10_dates$X131298.0.0 <- as.Date(icd10_dates$X131298.0.0)
icd10_dates$X131300.0.0 <- as.Date(icd10_dates$X131300.0.0)
icd10_dates$X131302.0.0 <- as.Date(icd10_dates$X131302.0.0)
icd10_dates$X131304.0.0 <- as.Date(icd10_dates$X131304.0.0)
icd10_dates$X131306.0.0 <- as.Date(icd10_dates$X131306.0.0)


icd10_dates$earliest_cad_date <- pmin(icd10_dates$X131298.0.0, icd10_dates$X131300.0.0,
                                      icd10_dates$X131302.0.0, icd10_dates$X131304.0.0,
                                      icd10_dates$X131306.0.0, na.rm = TRUE)

sum(is.na(icd10_dates$earliest_cad_date))

names(icd10_dates)[names(icd10_dates) == "eid"] <- "IID"







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


################################################################################
##### all well and good having the phenotypes but we need dates for all ########
##### of these ICD10 and ICD9 codes                                     ########
################################################################################


phenotypes_all <- merge(phenotypes_collapsed, phenotypes2, by = "eid")



names(phenotypes_all)[names(phenotypes_all) == "X30850.0.0"] <- "T"
names(phenotypes_all)[names(phenotypes_all) == "X41272.0.0"] <- "OPS"
names(phenotypes_all)[names(phenotypes_all) == "X20002.0.0"] <- "Self-report"
names(phenotypes_all)[names(phenotypes_all) == "eid"] <- "IID"



phenotypes_all <- merge(phenotypes_all, icd10_dates, by = "IID", all.x =TRUE)

phenotypes_all <- subset(phenotypes_all, select = -c(X131306.0.0, X131304.0.0, 
                                                     X131302.0.0, X131300.0.0, 
                                                     X131298.0.0))

####### LOCATING AND DATING OPERATIONS #########################################


phenotypes_all <- phenotypes_all %>%
  mutate("CAD_OP" = if_else(grepl("K40|K401|K402|K403|K404|K408|K409|K41|K411|
  K412|K413|K414|K418|K419|K42|K421|K422|K423|K424|K428|K429|K43|K431|K432|K433|
  K434|K438|K439|K44|K441|K442|K448|K449|K45|K451|K452|K453|K454|K455|K456|K458|
  K459|K46|K461|K462|K463|K464|K465|K468|K469|K49|K491|K492|K49.3|K494|K498|K499|
  K50|K501|K502|K50.3|K504|K508|K509|K75|K75.1|K75.2|K75.3|K75.4|K75.8|K75.9", phenotypes_all$OPS), 1, 0))


op_dates <- read.csv("male_ops_full.csv")


patterns <- c("K40", "K401", "K402", "K403", "K404", "K408", "K409", "K41", "K411",
              "K412", "K413", "K414", "K418", "K419", "K42", "K421", "K422", "K423", 
              "K424", "K428", "K429", "K43", "K431", "K432", "K433", "K434", "K438", 
              "K439", "K44", "K441", "K442", "K448", "K449", "K45", "K451", "K452", 
              "K453", "K454", "K455", "K456", "K458", "K459", "K46", "K461", "K462", 
              "K463", "K464", "K465", "K468", "K469", "K49", "K491", "K492", "K493", 
              "K494", "K498", "K499", "K50", "K501", "K502", "K50.3", "K504", "K508", 
              "K509", "K75", "K75.1", "K75.2", "K75.4", "K75.8", "K75.9")



relevant_ops <- op_dates %>%
  filter(grepl(paste(patterns, collapse = "|"), oper4))


names(relevant_ops)[names(relevant_ops) == "eid"] <- "IID"

phenotypes_all <- merge(phenotypes_all, relevant_ops, by = "IID", all.x = TRUE)


phenotypes_all <- subset(phenotypes_all, select = -c(dnx_hesin_oper_id))





################# GATHERING ALL CASES OF CAD INTO A SINGLE BINARY PHENOTYPE ####

phenotypes_all$CAD_ICD10 <- as.numeric(phenotypes_all$CAD_ICD10)
phenotypes_all$CAD_ICD9 <- as.numeric(phenotypes_all$CAD_ICD9)
phenotypes_all$CAD_OP <- as.numeric(phenotypes_all$CAD_OP)


phenotypes_all$CADBIN <- as.numeric(rowSums(phenotypes_all[, c("CAD_ICD10", "CAD_OP", "CAD_ICD9")]) > 0)



phenotypes_all$earliest_cad_date_all <- pmin(phenotypes_all$earliest_cad_date,
                                             phenotypes_all$opdate, na.rm = TRUE)




sum(any(phenotypes_all$CADBIN == 0 & (phenotypes_all$earliest_cad_date_all)))
sum(any(phenotypes_all$CADBIN == 0 & !is.na(phenotypes_all$earliest_cad_date_all)))


surv <- merge(phenotypes_all, surv, by = "IID", all.x = TRUE)













########### CALCULATING SURVIVOR VARIABLES - DATES ETC. #########################


surv$LTFBIN <- ifelse(surv$LTF == "" ,0,1)

# adding censoring date
surv$censdate <- Sys.Date()
surv$censdate[surv$LTFBIN == 1] <- surv$LTF[surv$LTFBIN == 1]

# NEED TO GATHER ALL THE DIFFERENT TYPES OF CAD AND THEIR DATES 
# INTO A SINGLE CAD DATE COLUMN 


# Extract month and day
surv$cadmonth <- month(surv$earliest_cad_date_all)
surv$cadday <- day(surv$earliest_cad_date_all)
surv$cadyear <- year(surv$earliest_cad_date_all)
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
surv$testosterone_deficiency <- ifelse(surv$T.x < 12, 1, 0)

# need complete cases for testosterone deficiency 
surv <- surv[complete.cases(surv$T.x), ]


write.csv(surv, "CAD_SURV.csv", row.names = TRUE)

# NEED TO COME BACK TO THIS AND CHECK IF WE CAN INCORPORATE THE ICD-9 CODES
# BUT I DON'T THINK IT WILL MAKE MUCH DIFFERENCE BECAUSE THE EARLIEST DATE 
# OF CAD IS RECORDED REGARDLESS


#################################################################################
####################### COPY THIS SECTION INTO THE RAP ##########################
########################################################################################

library(tidyverse)
library(survival)



surv <- read.csv("CAD_SURV.csv")
surv <- surv %>% select(-c(1))
surv <- surv %>% select(-c(9))
surv <- surv %>% select(-c(9))
surv <- surv %>% select(-c(15))



names(surv)[names(surv) == "diag_icd10.x"] <- "All_ICD10_diags"
names(surv)[names(surv) == "diag_icd9"] <- "All_ICD9_diags"
names(surv)[names(surv) == "T.x"] <- "T"
names(surv)[names(surv) == "OPS"] <- "All_OPS"

qrisk1 <- read.csv("qriskpt1.csv")
qrisk2 <- read.csv("qriskpt2.csv")

names(qrisk1)[names(qrisk1) == "Participant.ID"] <- "IID"
names(qrisk2)[names(qrisk2) == "Participant.ID"] <- "IID"

qrisksurv <- merge(merge(qrisk1, qrisk2, by = "IID", all = TRUE), surv, by = "IID", all = TRUE)
names(qrisksurv)[c(28, 29)] <- c("BMI", "Date_SLE")
names(qrisksurv)[c(30, 31)] <- c("Date_MIGRAINE", "Date_RHEUMARTH")
names(qrisksurv)[c(32, 33)] <- c("Date_OTHARTH", "Date_JUVARTH")
names(qrisksurv)[c(34, 35)] <- c("Date_SCHIZ", "Date_BIP")
names(qrisksurv)[c(36, 37)] <- c("Date_RECURDEP", "Date_DEPEP")
names(qrisksurv)[c(38, 39)] <- c("Medication", "Non_canc_illness")
names(qrisksurv)[c(2, 3)] <- c("AGERECRUIT1", "DEP_INDEX")
names(qrisksurv)[c(4, 5)] <- c("ETHNICITY", "CHOLESTEROL")
names(qrisksurv)[c(6, 8)] <- c("CHOLESTEROL_MEASURE", "HDL")
names(qrisksurv)[c(9, 10)] <- c("SBP", "DIABETES")
names(qrisksurv)[c(11, 12)] <- c("DATE_TYPE1DIAB", "DATE_TYPE2DIAB")
names(qrisksurv)[c(13, 14)] <- c("EVER_SMOKED", "CURRENT_SMOKER")
names(qrisksurv)[c(15, 16)] <- c("NUM_CIGS_DAILY", "CHOLESTEROL_MED")

qrisksurv <- qrisksurv %>% select(-c(17))
qrisksurv <- qrisksurv %>% select(-c(7))
qrisksurv <- qrisksurv %>% select(-c(17))

names(qrisksurv)[c(16, 17)] <- c("MICROALB_IN_URINE", "CREATININE_IN_URINE")
names(qrisksurv)[c(18, 19)] <- c("FATH_ILL", "FATH_ALIVE")
names(qrisksurv)[c(20, 21)] <- c("MOTH_ILL", "SIBS_ILL")
names(qrisksurv)[c(22, 23)] <- c("NON_CANC_ILLNESS", "DATE_AFIB")
qrisksurv <- qrisksurv %>% select(-c(24))


# removing people with no info on age or entrance to study 
qrisksurv <- qrisksurv[!is.na(qrisksurv$AGERECRUIT), ]


# MAKING NEW VARIABLES TO SUMMARISE THE DIFFERENT COVARIATES 



# diabetes 
qrisksurv$TYPE1DBIN <- ifelse(qrisksurv$DATE_TYPE1DIAB == "" ,0,1)
qrisksurv$TYPE2DBIN <- ifelse(qrisksurv$DATE_TYPE2DIAB == "" ,0,1)
qrisksurv$TYPE1DBIN <- factor(qrisksurv$TYPE1DBIN, levels=c(0,1))
qrisksurv$TYPE2DBIN <- factor(qrisksurv$TYPE2DBIN, levels=c(0,1))

qrisksurv <- qrisksurv[, -which(names(qrisksurv) %in% c("DATE_TYPE1DIAB", "DATE_TYPE2DIAB"))]


# cholesterol medication 
qrisksurv$CholesterolBIN <- ifelse(qrisksurv$CHOLESTEROL_MED == "Cholesterol loweing medication" | 
                                     qrisksurv$CHOLESTEROL_MED =="Cholesterol lowering medication|Blood pressure medication|Insulin" |
                                     qrisksurv$CHOLESTEROL_MED == "Cholesterol lowering medication|Blood pressure medication" | 
                                     qrisksurv$CHOLESTEROL_MED == "Cholesterol lowering medication|Insulin" ,1,0)

qrisksurv$CholesterolBIN <- factor(qrisksurv$CholesterolBIN, levels=c(0,1))
# atrial fibrilation
qrisksurv$AFIBBIN1 <- ifelse(qrisksurv$NON_CANC_ILLNESS == "atrial fibrillation" ,1,0)
qrisksurv$AFIBBIN1 <- factor(qrisksurv$AFIBBIN1, levels=c(0,1))

# systemic lupus erythematosis 
qrisksurv$SLEBIN <- ifelse(qrisksurv$NON_CANC_ILLNESS == "systemic lupus erythematosis/sle" ,1,0)
qrisksurv$SLEBIN <- factor(qrisksurv$SLEBIN, levels=c(0,1))

# migraine
qrisksurv$MIGRAINEBIN1 <- ifelse(qrisksurv$NON_CANC_ILLNESS == "migraine" ,1,0)
qrisksurv$MIGRAINEBIN1 <- factor(qrisksurv$MIGRAINEBIN1, levels=c(0,1))

# mental illnesses 
qrisksurv$MENTALBIN <- ifelse(qrisksurv$NON_CANC_ILLNESS %in% c("depression", "deliberate self-harm/suicide attempt", "schizophrenia"), 1, 0)
qrisksurv$MENTALBIN <- factor(qrisksurv$MENTALBIN, levels=c(0,1))

# erectile dysfunction 
qrisksurv$EDBIN <- ifelse(qrisksurv$NON_CANC_ILLNESS == "erectile dysfunction / impotence" ,1,0)
qrisksurv$EDBIN <- factor(qrisksurv$EDBIN, levels=c(0,1))

table(qrisksurv$NON_CANC_ILLNESS=="erectile dysfunction / impotence")
table(qrisksurv$EDBIN)

# atypical antipsychotic medication 
qrisksurv$ANTIPSYCHOTICMEDBIN <- ifelse(qrisksurv$Medication %in% c("risperidone", 
                                                                    "quetiapine", "olanzapine", "aripiprazole", 
                                                                    "clozapine", "haloperidol", "amisulpride", 
                                                                    "lurasidone", "paliperidone", "netiquette", 
                                                                    "sertindole", "zotepine") ,1,0)
qrisksurv$ANTIPSYCHOTICMEDBIN <- factor(qrisksurv$ANTIPSYCHOTICMEDBIN, levels=c(0,1))

# steroids  
qrisksurv$STEROIDSBIN <- ifelse(qrisksurv$Medication %in% c("prednisolone", "betamethasone",
                                                            "dexamethasone", "hydrocortisone", 
                                                            "cortisone", "depo-medrone", "deflazacort", 
                                                            "efcortesol", "methylprednisolone",
                                                            "triamcinolone") ,1,0)

qrisksurv$STEROIDSBIN <- factor(qrisksurv$STEROIDSBIN, levels=c(0,1))


# arthritis 
qrisksurv$ARTHBIN <- ifelse(qrisksurv$Date_RHEUMARTH == "" ,0,1)
qrisksurv$ARTHBIN <- factor(qrisksurv$ARTHBIN, levels=c(0,1))


# bp treatment

qrisksurv$BPTREATBIN <- ifelse(qrisksurv$Medication %in% c("enalapril", "lisinopril",
                                                           "perindopril", "ramipril", 
                                                           "losartan", "valsartan", "olmesartan", 
                                                           "amlodipine", "felodipine", "nifedipine", "diltiazem", 
                                                           "verapami", "indapamide", "bendroflumethiazide", "atenolol", 
                                                           "bisoprolol") ,1,0)

qrisksurv$BPTREATBIN <- factor(qrisksurv$STEROIDSBIN, levels=c(0,1))



# ethnicity 


qrisksurv$ETHNICITY <- ifelse(qrisksurv$ETHNICITY %in% c("Do not know", "Prefer not to answer", "Other ethnic group"), "Unknown",
                              ifelse(qrisksurv$ETHNICITY %in% c("African", "Black or Black British"), "African",
                                     ifelse(qrisksurv$ETHNICITY %in% c("Any other Asian background", "Asian or Asian British", "Bangladeshi", "Chinese", "Indian", "Pakistani"), "Asian",
                                            ifelse(qrisksurv$ETHNICITY %in% c("Any other mixed background", "Mixed", "White and Black African", "White and Black Caribbean", "White and Asian"), "Mixed", 
                                                   ifelse(qrisksurv$ETHNICITY %in% c("White", "Any other white background"), "White",
                                                          ifelse(qrisksurv$ETHNICITY %in% c("British"), "British", 
                                                                 ifelse(qrisksurv$ETHNICITY %in% c("Irish"), "Irish", qrisksurv$ETHNICITY)
                                                          )
                                                   )
                                            )
                                     )
                              )
)


qrisksurv$ETHNICITY <- ifelse(qrisksurv$ETHNICITY %in% c("Any other Black background", "Caribbean", "African"), "Black", qrisksurv$ETHNICITY)

qrisksurv$ETHNICITY <- as.factor(qrisksurv$ETHNICITY)
table(qrisksurv$ETHNICITY)
qrisksurv$ETHNICITY <- relevel(qrisksurv$ETHNICITY, ref = "White")


# chronic kidney disease - using microalbumin levels 

qrisksurv$CKDBIN <- ifelse(qrisksurv$MICROALB_IN_URINE < 30 ,0,1)
qrisksurv$CKDBIN  <- factor(qrisksurv$CKDBIN, levels=c(0,1))


# smoking - this one is more tricky 

qrisksurv$exsmoker <- ifelse(qrisksurv$EVER_SMOKED == "Yes" & qrisksurv$CURRENT_SMOKER=="No" ,"2",NA)
qrisksurv$nonsmoker <- ifelse(qrisksurv$EVER_SMOKED == "No" ,"1",NA)

qrisksurv$NUM_CIGS_DAILY <- as.numeric(qrisksurv$NUM_CIGS_DAILY)

qrisksurv$SmokingCategory <- ifelse(qrisksurv$NUM_CIGS_DAILY < 10 & qrisksurv$NUM_CIGS_DAILY > 0, "3",
                                    ifelse(qrisksurv$NUM_CIGS_DAILY >= 10 & qrisksurv$NUM_CIGS_DAILY < 20, "4", 
                                           ifelse(qrisksurv$NUM_CIGS_DAILY >= 20, "5", NA)))



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



# family history 

qrisksurv$HDFATHBIN <- ifelse(grepl("Heart disease", qrisksurv$FATH_ILL), 1, "")
qrisksurv$HDMOTHBIN <- ifelse(grepl("Heart disease", qrisksurv$MOTH_ILL), 1, "")
qrisksurv$HDSIBBIN <- ifelse(grepl("Heart disease", qrisksurv$SIBS_ILL), 1, "")

qrisksurv$FAMHISTBIN <- paste0(qrisksurv$HDFATHBIN, qrisksurv$HDMOTHBIN, qrisksurv$HDSIBBIN)
qrisksurv$FAMHISTBIN <- ifelse(qrisksurv$FAMHISTBIN>0, 1, 0)

qrisksurv$FAMHISTBIN <- factor(qrisksurv$FAMHISTBIN, levels = c(0,1))



# deprivation

# Check for CAD incidence before entry into the study
qrisksurv$TIMETOCAD <- qrisksurv$timetoCAD - qrisksurv$AGERECRUIT
table(qrisksurv$TIMETOCAD)
qrisksurv <- qrisksurv[is.na(qrisksurv$TIMETOCAD) | qrisksurv$TIMETOCAD >= 0, ]


# gathering the relevant columns



qrisksurv <- select(qrisksurv, AGERECRUIT1, CADBIN, timetoEVENT, testosterone_deficiency, 
                    BMI, TYPE1DBIN, TYPE2DBIN, EDBIN, SLEBIN, AFIBBIN1, UKBBSMOKING, 
                    MENTALBIN, MIGRAINEBIN1, STEROIDSBIN, ANTIPSYCHOTICMEDBIN, 
                    CholesterolBIN, ARTHBIN, BPTREATBIN, CKDBIN, CHOLESTEROL_MEASURE, 
                    HDL, SBP, FAMHISTBIN, ETHNICITY, DEP_INDEX)


qrisksurv <- qrisksurv[complete.cases(qrisksurv), ]


# going to need to do imputation 


### ADJUSTING FOR AGE AND ERECTILE DYSFUNCTION 

CADsurv <- Surv(qrisksurv$timetoEVENT, qrisksurv$CADBIN)


qrisksurv$testosterone_deficiency <- factor(qrisksurv$testosterone_deficiency, levels=c(0,1))
km1 <- survfit(CADsurv~qrisksurv$testosterone_deficiency)


summary(km1, censored = T)

plot(km1, lty = c(1, 2), col = c("blue", "red"),
     xlab = "Time", ylab = "Survival Probability",
     main = "Kaplan-Meier Curves for Testosterone Deficiency")
legend("bottomright", legend = c("Non-Testosterone Deficient", "Testosterone Deficient"),
       col = c("blue", "red"), lty = c(1, 2))




coxED <- coxph(CADsurv~qrisksurv$testosterone_deficiency)
summary(coxED)







# ADJUSTING FOR AGE, ED AND DIABETES

CADsurv <- Surv(qrisksurv$timetoEVENT, qrisksurv$CADBIN)


coxEDD <- coxph(CADsurv~qrisksurv$testosterone_deficiency+qrisksurv$EDBIN)
summary(coxEDD)

sum(qrisksurv$EDBIN==1)




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




coxALL <- coxph(CADsurv~qrisksurv$TESTOSTERONE_DEFICIENCY+qrisksurv$AGERECRUIT+qrisksurv$EDBIN+
                  qrisksurv$TYPE1DBIN + qrisksurv$TYPE2DBIN + qrisksurv$CholesterolBIN + qrisksurv$AFIBBIN1
                + qrisksurv$SLEBIN + qrisksurv$MIGRAINEBIN1 +qrisksurv$MENTALBIN +
                  qrisksurv$UKBBSMOKING + qrisksurv$STEROIDSBIN + qrisksurv$ANTIPSYCHOTICMEDBIN +
                  qrisksurv$BMI1 + qrisksurv$ARTHBIN + qrisksurv$BPTREATBIN + qrisksurv$CKDBIN +
                  qrisksurv$Cholesterol_level + qrisksurv$HDL_instance0 + qrisksurv$SBP + 
                  qrisksurv$FAMHISTBIN +qrisksurv$Ethnicity +qrisksurv$Deprivation, data = qrisksurv)



summary(coxALL)
coxALL


