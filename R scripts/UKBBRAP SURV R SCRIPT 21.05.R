

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


# reading in the file which has all the covariates that we want to control for 
# including the survivor variables like age of recruitment etc. 

surv <- read.csv("data_participant_surv_2.csv")

colnames(surv) <- c("IID", "T", "CAD", "AGERECRUIT", "MONTHBIRTH", "YEARBIRTH", "LTF", "DATEASSESSMENT")



################################################################################
######################### DATING CAD PHENOTYPES  ###############################
################################################################################

## read in data on ICD_10 codes 
## icd10_1 is the participants and the dates of their diagnosis of specific 
## types of conditions 
## this dates them by category of condition rather than the very specific 
## condition types that are listed in the ICD_10 codes themselves 
## icd10_2 is a row per diagnosis of each individual 
## so each individual may have as many rows as they have conditions
## these are described by their ICD_10 code


icd10_1 <- read.csv("male_icd10_1.csv")
icd10_2 <- read.csv("male_icd10_2.csv")


# assign the patterns of characters that we want to search for in the 
# ICD_10 codes 
# all of the cardiovascular conditions begin with an I and the 
# subsequent letters describe the more specific subtypes 


patterns <- c("I21", "I210", "I211", "I212", "I213", "I214", "I219", "I21X",
              "I22", "I220", "I221", "I228", "I229", "I23", "I23.1", "I23.2",
              "I23.3", "I23.4", "I23.5", "I23.6", "I238", "I24", "I240",
              "I241", "I248", "I249", "I252")


# now we are filtering the long list of conditions which has a row 
# for each individual and just keeping the ones with the cardiovascular 
# ICD_10 patterns 


relevant_icd10 <- icd10_2 %>%
  filter(grepl(paste(patterns, collapse = "|"), diag_icd10))


# then merge with the file which has all the dates to make it possible 
# to date the cardiovascular conditions 
# the different columns from the ICD_10_1 file represent 
# the different types of conditions 
# within them is the date which the person was diagnosed with the condition 


icd10_dates <- merge(icd10_1, relevant_icd10, by = "eid")
icd10_dates$X131298.0.0 <- as.Date(icd10_dates$X131298.0.0)
icd10_dates$X131300.0.0 <- as.Date(icd10_dates$X131300.0.0)
icd10_dates$X131302.0.0 <- as.Date(icd10_dates$X131302.0.0)
icd10_dates$X131304.0.0 <- as.Date(icd10_dates$X131304.0.0)
icd10_dates$X131306.0.0 <- as.Date(icd10_dates$X131306.0.0)


# now it is possible to find the minimum value of the dates of the 
# different diagnoses of types of CAD conditions 
# so we can use the pmin function to select the earliest date
# this date is required as we want the earliest diagnosis of CAD 
# to be used as the date of the event in the survival analysis 


icd10_dates$earliest_cad_date <- pmin(icd10_dates$X131298.0.0, icd10_dates$X131300.0.0,
                                      icd10_dates$X131302.0.0, icd10_dates$X131304.0.0,
                                      icd10_dates$X131306.0.0, na.rm = TRUE)


# checking if there are any NAs to see if anything has gone wrong in 
# finding the minimum value 
# now for every CAD condition recorded there should be the date associated with 
# the earliest diagnosis 

sum(is.na(icd10_dates$earliest_cad_date))

# changing the name of the ID column so it is possible to merge with the next 
# table

names(icd10_dates)[names(icd10_dates) == "eid"] <- "IID"







################################################################################
############# phenotyping CAD cases ############################################
################################################################################


## phenotypes 1 contains the list of participants and their icd10 and icd9
## diagnoses 
## there are very few ICD9 diagnoses because they are an older form
## i believe ICD9 codes were just used in scotland
## phenotypes 2 has testosterone levels 
## also has all of the recorded operations for these individuals 
## and any self reported illness 

phenotypes <- read.csv("new_cad_1.csv")
phenotypes2 <- read.csv("new_cad_2.csv")


# here we are going to collapse the icd10 data so that there is 
# not one row per condition per individual but so they are all 
# in one row for that individual

phenotypes_collapsed <- phenotypes %>%
  group_by(eid) %>%
  summarize(
    diag_icd10 = paste(diag_icd10, collapse = ", "),
    diag_icd9 = paste(diag_icd9, collapse = "")
  ) %>%
  ungroup()


# then we are going to pull out any CAD related conditions based on their 
# icd10 code and create a new column called CAD_ICD10 which places a 1
# if any CAD conditions were present in their ICD10 list 
# and a 0 if there were not any 

phenotypes_collapsed <- phenotypes_collapsed %>%
  mutate(CAD_ICD10 = if_else(grepl("I21|I210|I211|I212|I213|
                            I214|I219|I21X|I22|I220|I221|I228|I229|
                            I23|I23.1|I23.2|I23.3|I23.4|I23.5|
                            I23.6|I238|I24|I240|I241|I248|I249|
                            I252", phenotypes_collapsed$diag_icd10), 1, 0))


# now doing the same for ICD9 codes 
# the ones listed in this code are the way CAD is recorded in ICD9


phenotypes_collapsed <- phenotypes_collapsed %>%
  mutate("CAD_ICD9" = if_else(grepl("410|411|
                                    412|413|414", phenotypes_collapsed$diag_icd9), 1, 0))


## now we have coded whether or not individuals have or do not have CAD 
## as defined by a long list of ICD10 and ICD9 codes 
## now we are going to merge this with the file which has data on operations 
## and self reported CAD 


phenotypes_all <- merge(phenotypes_collapsed, phenotypes2, by = "eid")

## and rename the columns so we can understand them 

names(phenotypes_all)[names(phenotypes_all) == "X30850.0.0"] <- "T"
names(phenotypes_all)[names(phenotypes_all) == "X41272.0.0"] <- "OPS"
names(phenotypes_all)[names(phenotypes_all) == "X20002.0.0"] <- "Self-report"
names(phenotypes_all)[names(phenotypes_all) == "eid"] <- "IID"



## adding the ICD_10 dates file onto the phenotypes_all file 

phenotypes_all <- merge(phenotypes_all, icd10_dates, by = "IID", all.x =TRUE)


## removing all the columns for the different types of CAD 
## as we do not need these anymore


phenotypes_all <- subset(phenotypes_all, select = -c(X131306.0.0, X131304.0.0, 
                                                     X131302.0.0, X131300.0.0, 
                                                     X131298.0.0))

####### LOCATING AND DATING OPERATIONS #########################################


## now we are going to pull out all of the operations which are associated 
## with CAD - these fall under this long list of codes 
## if people have had these operations, they get a 1, if not they get a 0


phenotypes_all <- phenotypes_all %>%
  mutate("CAD_OP" = if_else(grepl("K40|K401|K402|K403|K404|K408|K409|K41|K411|
  K412|K413|K414|K418|K419|K42|K421|K422|K423|K424|K428|K429|K43|K431|K432|K433|
  K434|K438|K439|K44|K441|K442|K448|K449|K45|K451|K452|K453|K454|K455|K456|K458|
  K459|K46|K461|K462|K463|K464|K465|K468|K469|K49|K491|K492|K49.3|K494|K498|K499|
  K50|K501|K502|K50.3|K504|K508|K509|K75|K75.1|K75.2|K75.3|K75.4|K75.8|K75.9", phenotypes_all$OPS), 1, 0))


## now we are going to find out what the dates of these operations were
## this is using a similar method as we used for the ICD10 data


## this file has data on the operations of each individual and when these operations
## happened

op_dates <- read.csv("male_ops_full.csv")


patterns <- c("K40", "K401", "K402", "K403", "K404", "K408", "K409", "K41", "K411",
              "K412", "K413", "K414", "K418", "K419", "K42", "K421", "K422", "K423", 
              "K424", "K428", "K429", "K43", "K431", "K432", "K433", "K434", "K438", 
              "K439", "K44", "K441", "K442", "K448", "K449", "K45", "K451", "K452", 
              "K453", "K454", "K455", "K456", "K458", "K459", "K46", "K461", "K462", 
              "K463", "K464", "K465", "K468", "K469", "K49", "K491", "K492", "K493", 
              "K494", "K498", "K499", "K50", "K501", "K502", "K50.3", "K504", "K508", 
              "K509", "K75", "K75.1", "K75.2", "K75.4", "K75.8", "K75.9")


## now filtering the operations data to keep only the operations associated 
## with CAD


relevant_ops <- op_dates %>%
  filter(grepl(paste(patterns, collapse = "|"), oper4))


## renaming the ID file so we can merge with the phenotypes file

names(relevant_ops)[names(relevant_ops) == "eid"] <- "IID"

## merging the operations dates to the phenotypes_all file

phenotypes_all <- merge(phenotypes_all, relevant_ops, by = "IID", all.x = TRUE)

## removing redundant columns 

phenotypes_all <- subset(phenotypes_all, select = -c(dnx_hesin_oper_id))





########### GATHERING ALL CASES OF CAD INTO A SINGLE BINARY PHENOTYPE ##########


## setting all the CAD binary outcomes as numerics 

phenotypes_all$CAD_ICD10 <- as.numeric(phenotypes_all$CAD_ICD10)
phenotypes_all$CAD_ICD9 <- as.numeric(phenotypes_all$CAD_ICD9)
phenotypes_all$CAD_OP <- as.numeric(phenotypes_all$CAD_OP)


## telling R that if there is a 1 in any of these 3 columns to put a
## 1 in our new CADBIN column 
## this CADBIN column has a 1 if CAD has been identified by either ICD10, ICD9
## or operations codes 


phenotypes_all$CADBIN <- as.numeric(rowSums(phenotypes_all[, c("CAD_ICD10", "CAD_OP", "CAD_ICD9")]) > 0)




## now we are merging the earliest CAD date column from the ICD_10 data 
## and the operation date for those who had the operation, and selecting the 
## first instance

phenotypes_all$earliest_cad_date_all <- pmin(phenotypes_all$earliest_cad_date,
                                             phenotypes_all$opdate, na.rm = TRUE)



## then checking if there is anyone that does not have cad and has a date suggesting 
## they have cad 

sum(any(phenotypes_all$CADBIN == 0 & !is.na(phenotypes_all$earliest_cad_date_all)))

## now checking if there is anyone that does have cad but does not have a date 

sum(any(phenotypes_all$CADBIN == 1 & is.na(phenotypes_all$earliest_cad_date_all)))




## now merging our CAD survivorship info with the surv file which we 
## loaded in first and contains all the relevant covariates 


surv <- merge(phenotypes_all, surv, by = "IID", all.x = TRUE)


########### CALCULATING SURVIVOR VARIABLES - DATES ETC. #########################

# if they have a date in their lost to follow up column, place a 1
# otherwise, leave as 0 - now we have a binary column which says whether 
# someone was lost to follow up or not

surv$LTFBIN <- ifelse(surv$LTF == "" ,0,1)

# adding censoring date as the current date 
surv$censdate <- Sys.Date()
surv$censdate[surv$LTFBIN == 1] <- surv$LTF[surv$LTFBIN == 1]

# Extract month and day of the earliest recorded CAD instance 
surv$cadmonth <- month(surv$earliest_cad_date_all)
surv$cadday <- day(surv$earliest_cad_date_all)
surv$cadyear <- year(surv$earliest_cad_date_all)
surv$censyear <- year(surv$censdate)

# changing months to numbers for month of birth
surv$monthbirthnum <- as.integer(factor(surv$MONTHBIRTH, levels = month.name))

# creating a censoring variable = anyone who does not come up as a CAD case
surv$censored <- ifelse(surv$CADBIN == 0, 1, 0)

# creating a year of recruitment variable - this will be inaccurate 
surv$year_of_recruitment <- surv$YEARBIRTH + surv$AGERECRUIT

# creating a time to event variable
surv$timetoCAD <- surv$cadyear-surv$YEARBIRTH

# creating a time to censoring variable 
surv$timetoCENSOR <- ifelse(surv$censored == 1, surv$censyear - surv$YEARBIRTH, NA)

table(surv$censdate)

# creating a general time to event variable 
surv$timetoEVENT <- ifelse(is.na(surv$timetoCENSOR), surv$timetoCAD, surv$timetoCENSOR)

# create a variable for testosterone deficiency 
surv$testosterone_deficiency <- ifelse(surv$T.x < 12, 1, 0)

# need complete cases for testosterone deficiency 
surv <- surv[complete.cases(surv$T.x), ]


# writing this surv file out so it can be used in the RAP to run the actual
# models

write.csv(surv, "CAD_SURV", row.names = TRUE)

# NEED TO COME BACK TO THIS AND CHECK IF WE CAN INCORPORATE THE ICD-9 CODES
# BUT I DON'T THINK IT WILL MAKE MUCH DIFFERENCE BECAUSE THE EARLIEST DATE 
# OF CAD IS RECORDED REGARDLESS


#################################################################################
####################### COPY THIS SECTION INTO THE RAP ##########################
########################################################################################

library(tidyverse)
library(survival)

# reading in the surv file that we have just written out 
# and removing some redundant columns 

surv <- read.csv("CAD_SURV.csv")
surv <- surv %>% select(-c(1))
surv <- surv %>% select(-c(9))
surv <- surv %>% select(-c(9))
surv <- surv %>% select(-c(15))

# renaming some of the columns

names(surv)[names(surv) == "diag_icd10.x"] <- "All_ICD10_diags"
names(surv)[names(surv) == "diag_icd9"] <- "All_ICD9_diags"
names(surv)[names(surv) == "T.x"] <- "T"
names(surv)[names(surv) == "OPS"] <- "All_OPS"

# reading in the covariates file

qrisk1 <- read.csv("qriskpt1.csv")
qrisk2 <- read.csv("qriskpt2.csv")

# renaming the ID column for merging

names(qrisk1)[names(qrisk1) == "Participant.ID"] <- "IID"
names(qrisk2)[names(qrisk2) == "Participant.ID"] <- "IID"

# merging and renaming the rest of the columns 

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

# removing some redundant columns 

qrisksurv <- qrisksurv %>% select(-c(17))
qrisksurv <- qrisksurv %>% select(-c(7))
qrisksurv <- qrisksurv %>% select(-c(17))

# renaming the rest of the columns

names(qrisksurv)[c(16, 17)] <- c("MICROALB_IN_URINE", "CREATININE_IN_URINE")
names(qrisksurv)[c(18, 19)] <- c("FATH_ILL", "FATH_ALIVE")
names(qrisksurv)[c(20, 21)] <- c("MOTH_ILL", "SIBS_ILL")
names(qrisksurv)[c(22, 23)] <- c("NON_CANC_ILLNESS", "DATE_AFIB")
qrisksurv <- qrisksurv %>% select(-c(24))


# removing people with no info on age or entrance to study 
qrisksurv <- qrisksurv[!is.na(qrisksurv$AGERECRUIT), ]


########  MAKING NEW VARIABLES TO SUMMARISE THE DIFFERENT COVARIATES ##########



# diabetes 
qrisksurv$TYPE1DBIN <- ifelse(qrisksurv$DATE_TYPE1DIAB == "" ,0,1)
qrisksurv$TYPE2DBIN <- ifelse(qrisksurv$DATE_TYPE2DIAB == "" ,0,1)
qrisksurv$TYPE1DBIN <- factor(qrisksurv$TYPE1DBIN, levels=c(0,1))
qrisksurv$TYPE2DBIN <- factor(qrisksurv$TYPE2DBIN, levels=c(0,1))

table(qrisksurv$TYPE1DBIN)
table(qrisksurv$TYPE2DBIN)

qrisksurv <- qrisksurv[, -which(names(qrisksurv) %in% c("DATE_TYPE1DIAB", "DATE_TYPE2DIAB"))]


# cholesterol medication 
qrisksurv$CholesterolBIN <- ifelse(qrisksurv$CHOLESTEROL_MED == "Cholesterol lowering medication" | 
                                     qrisksurv$CHOLESTEROL_MED =="Cholesterol lowering medication|Blood pressure medication|Insulin" |
                                     qrisksurv$CHOLESTEROL_MED == "Cholesterol lowering medication|Blood pressure medication" | 
                                     qrisksurv$CHOLESTEROL_MED == "Cholesterol lowering medication|Insulin" ,1,0)

qrisksurv$CholesterolBIN <- factor(qrisksurv$CholesterolBIN, levels=c(0,1))
table(qrisksurv$CHOLESTEROL_MED, qrisksurv$CholesterolBIN)


# atrial fibrilation
qrisksurv$AFIBBIN1 <- ifelse(qrisksurv$NON_CANC_ILLNESS == "atrial fibrillation" ,1,0)
qrisksurv$AFIBBIN1 <- factor(qrisksurv$AFIBBIN1, levels=c(0,1))

table(qrisksurv$AFIBBIN1, qrisksurv$NON_CANC_ILLNESS=="atrial fibrillation")

# systemic lupus erythematosis 
qrisksurv$SLEBIN <- ifelse(qrisksurv$NON_CANC_ILLNESS == "systemic lupus erythematosis/sle" ,1,0)
qrisksurv$SLEBIN <- factor(qrisksurv$SLEBIN, levels=c(0,1))

table(qrisksurv$SLEBIN, qrisksurv$NON_CANC_ILLNESS=="systemic lupus erythematosis/sle")
table(qrisksurv$NON_CANC_ILLNESS)

# need to check if there is a more detailed illness column 

# migraine
qrisksurv$MIGRAINEBIN1 <- ifelse(qrisksurv$NON_CANC_ILLNESS == "migraine" ,1,0)
qrisksurv$MIGRAINEBIN1 <- factor(qrisksurv$MIGRAINEBIN1, levels=c(0,1))
table(qrisksurv$MIGRAINEBIN1, qrisksurv$NON_CANC_ILLNESS=="migraine")

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



qrisksurv <- select(qrisksurv, AGERECRUIT1, CADBIN, T, timetoEVENT, testosterone_deficiency, 
                    BMI, TYPE1DBIN, TYPE2DBIN, EDBIN, SLEBIN, AFIBBIN1, UKBBSMOKING, 
                    MENTALBIN, MIGRAINEBIN1, STEROIDSBIN, ANTIPSYCHOTICMEDBIN, 
                    CholesterolBIN, ARTHBIN, BPTREATBIN, CKDBIN, CHOLESTEROL_MEASURE, 
                    HDL, SBP, FAMHISTBIN, ETHNICITY, DEP_INDEX)


qrisksurv <- qrisksurv[complete.cases(qrisksurv), ]


# CREATING AGE GROUP VARIABLES AND ASSIGNED T DEFICIENCIES 


qrisksurv$age_group <- cut(qrisksurv$AGERECRUIT1,
                           breaks = c(-Inf, 50, 60, 70, Inf),
                           labels = c("40-50", "50-60", "60-70", "70+"),
                           right = FALSE)



t_deficient_age <- function(age_group, T) {
  if (age_group == "40-50") {
    return(ifelse(T < 8.74, 1, 0))
  } else if (age_group == "50-60") {
    return(ifelse(T < 7.47, 1, 0))
  } else if (age_group == "60-70") {
    return(ifelse(T < 6.80, 1, 0))
  } else if (age_group == "70+") {
    return(ifelse(T < 5.41, 1, 0))
  } else {
    return(NA)
  }
}


qrisksurv$testosterone_deficient_age <- mapply(t_deficient_age, qrisksurv$age_group, qrisksurv$T)






# going to need to do imputation 


### ADJUSTING FOR AGE AND ERECTILE DYSFUNCTION 

CADsurv <- Surv(qrisksurv$timetoEVENT, qrisksurv$CADBIN)


qrisksurv$testosterone_deficient_age <- factor(qrisksurv$testosterone_deficient_age, levels=c(0,1))
km1 <- survfit(CADsurv~qrisksurv$testosterone_deficient_age)


summary(km1, censored = T)

plot(km1, lty = c(1, 2), col = c("blue", "red"),
     xlab = "Time", ylab = "Survival Probability",
     main = "Kaplan-Meier Curves for Testosterone Deficiency")
legend("bottomright", legend = c("Non-Testosterone Deficient", "Testosterone Deficient"),
       col = c("blue", "red"), lty = c(1, 2))




cox_model <- coxph(CADsurv~qrisksurv$testosterone_deficient_age)
summary(cox_model)



cox.zph_test <- cox.zph(cox_model)
print(cox.zph_test)
plot(cox.zph_test)

install.packages("survminer")
library(survminer)

# Plot diagnostic plots using survminer
ggcoxzph(cox.zph_test)


# ADJUSTING FOR ED 

CADsurv <- Surv(qrisksurv$timetoEVENT, qrisksurv$CADBIN)


coxEDD <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$EDBIN)
summary(coxEDD)

sum(qrisksurv$EDBIN==1)



# ADJUSTING FOR OTHER THINGS

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$BMI)
summary(cox) # ASSOCIATION ATTENUATES


cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$TYPE2DBIN)
summary(cox)

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$SLEBIN)
summary(cox)

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$AFIBBIN1)
summary(cox)

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$UKBBSMOKING)
summary(cox)

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$MENTALBIN)
summary(cox)

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$MIGRAINEBIN1)
summary(cox)

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$STEROIDSBIN)
summary(cox)

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$ANTIPSYCHOTICMEDBIN)
summary(cox)

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$CholesterolBIN)
summary(cox)

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$ARTHBIN)
summary(cox)

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$BPTREATBIN)
summary(cox)

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$CKDBIN)
summary(cox)

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$CHOLESTEROL_MEASURE)
summary(cox)

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$HDL)
summary(cox)

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$SBP)
summary(cox)

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$FAMHISTBIN)
summary(cox)

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$ETHNICITY)
summary(cox)

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$DEP_INDEX)
summary(cox)




# ADJUSTING FOR EVERYTHING

cox <- coxph(CADsurv~qrisksurv$testosterone_deficient_age+qrisksurv$DEP_INDEX +
               qrisksurv$BMI + qrisksurv$TYPE1DBIN + qrisksurv$TYPE2DBIN +
               qrisksurv$AFIBBIN1 + qrisksurv$UKBBSMOKING + 
               qrisksurv$MENTALBIN + qrisksurv$MIGRAINEBIN1 + qrisksurv$STEROIDSBIN +
               qrisksurv$ANTIPSYCHOTICMEDBIN + qrisksurv$CholesterolBIN + qrisksurv$ARTHBIN +
               qrisksurv$BPTREATBIN + qrisksurv$CKDBIN + qrisksurv$CHOLESTEROL_MEASURE + 
               qrisksurv$HDL + qrisksurv$SBP + qrisksurv$FAMHISTBIN + qrisksurv$ETHNICITY)
summary(cox)



















# having a look at the distribution of BMI AND TESTOSTERONE ETC. 


ggplot(qrisksurv, aes(x = BMI, fill = factor(testosterone_deficient_age))) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of BMI Values by Testosterone Deficiency Status",
       x = "BMI",
       fill = "Testosterone Deficiency Age") +
  scale_fill_manual(values = c("0" = "blue", "1" = "red"), labels = c("No Deficiency", "Deficiency")) +
  theme_minimal()


# looking at the relationship between BMI and testosterone in each 
# age group 


ggplot(qrisksurv, aes(x = T, y = BMI, color = age_group, group = age_group)) +
  geom_point() +
  labs(title = "Association between Testosterone Levels and BMI by Age Group",
       x = "Testosterone Levels (nmol/L)",
       y = "BMI",
       color = "Age Group") +
  theme_minimal()




results <- lapply(levels(qrisksurv$age_group), function(group) {
  subset_data <- subset(qrisksurv, age_group == group)
  model <- lm(BMI ~ T, data = subset_data)
  summary(model)
})



# Print results
names(results) <- levels(qrisksurv$age_group)
results





ggplot(qrisksurv, aes(x = T, y = BMI, color = age_group)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Association between Testosterone Levels and BMI by Age Group",
       x = "Testosterone Levels (nmol/L)",
       y = "BMI",
       color = "Age Group") +
  theme_minimal()



