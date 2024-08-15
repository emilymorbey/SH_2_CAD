library(tidyverse)
library(dplyr)
library(survival)
library(lubridate)
library(stringr)
library(cli)
library(data.table)
library(ggplot2)


#############            STEP 1: DATA CLEANING              #####################  

################################################################################
###########  THIS SECTION OF THE SCRIPT WAS RUN IN RSTUDIO  #####################
##########   ON A LOCAL DEVICE. DATA WAS DOWNLOADED USING  #####################
##########   COHORT BROWSER IN UK BIOBANK RAP. DATA ON     #####################
#########    SURVIVAL VARIABLES, DIAGNOSTIC CRITERIA FOR  ######################
##########   CAD AND DATA ON COVARIATES WERE DOWNLOADED    #####################
##########   INTO SEPARATE FILES ONTO A LOCAL DRIVE        #####################
###########  FILE NAMES FOR EACH OF THESE FILE TYPES       #####################
##########   ARE AS FOLLOWS:                                      ###############
##########   SURVIVAL DATA: data_participant_surv_2.csv          ################
##########   ICD10 CODES:   male_icd10_1.csv and male_icd10_2.csv   ############
#########    POTENTIAL CAD PHENOTYPES: new_cad_1.csv and new_cad_2.csv #########
#########    DATE OF ATTENDANCE: attendance_participant.csv            #########
#########    MEDICATIONS: medications_participant.csv                  #########
#########    SOCIODEMOGRAPHICS: sociodemographics_participant.csv      #########
#########    OTHER DISEASES: other_diseases_participant.csv            #########
#########    DIABETES: diabetes3_participant.csv                       #########
########                                                               #########
########     these variables were downloaded separately owing to       #########
########     limited capacity of UKBRAP to store lots of variables     #########
########     at once in a data table                                   #########
#################################################################################




setwd("C:/Users/emorb/OneDrive - University of Cambridge/PhD/MR/Testosterone_CAD_MR/Testosterone CAD MR R files")


################################################################################
##################### READING IN SURVIVAL COVARIATES  ##########################
################################################################################


# reading in the file which has survivor variables like age of recruitment etc. 

library(tidyverse)
library(dplyr)
library(survival)
library(lubridate)
library(stringr)
library(cli)
library(data.table)
library(ggplot2)

setwd("C:/Users/emorb/OneDrive - University of Cambridge/PhD/MR/Testosterone_CAD_MR/Testosterone CAD MR R files")


################################################################################
##################### ADDING IN COVARIATES #####################################
################################################################################


# reading in the file which has all the covariates that we want to control for 
# including the survivor variables like age of recruitment etc. 

surv <- read.csv("TestosteroneCAD/data_participant_surv_2.csv")

colnames(surv) <- c("IID", "T", "CAD", "AGERECRUIT", "MONTHBIRTH", "YEARBIRTH", "LTF", "DATEASSESSMENT")

surv <- surv %>% select(-CAD)



################################################################################
######################### DATING AND PHENOTYPING CAD PHENOTYPES  ###############
################################################################################

## read in data on ICD_10 codes 
## icd10_1 is the participants and the dates of their diagnosis of specific 
## types of conditions 
## this dates them by category of condition rather than the very specific 
## condition types that are listed in the ICD_10 codes themselves 
## icd10_2 is a row per diagnosis of each individual 
## so each individual may have as many rows as they have conditions
## these are described by their ICD_10 code


icd10_1 <- read.csv("TestosteroneCAD/male_icd10_1.csv")
icd10_2 <- read.csv("TestosteroneCAD/male_icd10_2.csv")


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

phenotypes <- read.csv("TestosteroneCAD/new_cad_1.csv")
phenotypes2 <- read.csv("TestosteroneCAD/new_cad_2.csv")


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

op_dates <- read.csv("TestosteroneCAD/male_ops_full.csv")


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

phenotypes_condensed <- phenotypes_all %>% select("IID", "T", "CADBIN", "earliest_cad_date_all")


## now merging our CAD survivorship info with the surv file which we 
## loaded in first and contains all the relevant covariates 


surv <- merge(phenotypes_condensed, surv, by = "IID", all.x = TRUE)
names(surv)[names(surv) == "T.x"] <- "T"
surv <- surv %>% select(-"T.y")

colnames(surv)


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

surv_with_dates <- surv

# changing months to numbers for month of birth
surv$monthbirthnum <- as.integer(factor(surv$MONTHBIRTH, levels = month.name))

# creating a censoring variable = anyone who does not come up as a CAD case
surv$censored <- ifelse(surv$CADBIN == 0, 1, 0)


## DATE OF ASSESSMENT VARIABLE

dates <- read.csv("TestosteroneCAD/attendance_participant.csv")

names(dates)[names(dates) == "Participant.ID"] <- "IID"

surv <- merge(surv, dates, by = "IID")

names(surv)[names(surv) == "Date.of.attending.assessment.centre...Instance.0"] <- "ASSESSMENT_DATE"


# creating a year of recruitment variable - this will be inaccurate 
surv$year_of_recruitment <- year(surv$ASSESSMENT_DATE)

# creating a time to event variable
surv$timetoCAD <- surv$cadyear-surv$year_of_recruitment

# creating a time to censoring variable 
surv$timetoCENSOR <- ifelse(surv$censored == 1, surv$censyear - surv$year_of_recruitment, NA)

# creating a general time to event variable 
surv$timetoEVENT <- ifelse(is.na(surv$timetoCENSOR), surv$timetoCAD, surv$timetoCENSOR)

# create a variable for testosterone deficiency 
# surv$testosterone_deficiency <- ifelse(surv$T.x < 12, 1, 0)

# need complete cases for testosterone deficiency 
# surv <- surv[complete.cases(surv$T.x), ]

any(duplicated(surv$IID))
any(duplicated(surv$IID))
surv <- surv[!duplicated(surv$IID), ]
surv <- surv[!duplicated(surv$IID), ]

surv$timetoCAD2 <- difftime(surv$earliest_cad_date_all, surv$timetoEVENT, units = "weeks")


# select relevant columns 

surv_key_variables <- surv %>% select("IID", "T", "CADBIN",
                                      "timetoEVENT")


# writing this surv file out so it can be used in the RAP to run the actual
# models

write.csv(surv_key_variables, "CAD_SURV.csv", row.names = TRUE)




library(tidyverse)
library(survival)

# reading in the surv file that we have just written out 
# and removing some redundant columns 

surv <- read.csv("CAD_SURV.csv")



################################################################################
################# ADDING IN COVARIATES #########################################
################################################################################



###### ALL COVARIATES READ IN  ################################################

medications <- read.csv("TestosteroneCAD/medications_participant.csv")
sociodemographics <- read.csv("TestosteroneCAD/sociodemographics_participant.csv")
other_illness <- read.csv("TestosteroneCAD/other_diseases_participant.csv")




###### DIABETES ################################################################

diabetes <- read.csv("TestosteroneCAD/diabetes3_participant.csv")


colnames(diabetes) <- c("IID", "SELFREPORT", "MEDICATION", "DOCTOR", "HBA1C", "ICD10", "ICD9")


#### type 1

diabetes <- diabetes %>%
  mutate(TYPE1DIAB = if_else(grepl("E10|O240", diabetes$ICD10) | 
                               grepl("1222", diabetes$SELFREPORT)|
                               grepl("25001|25011|25021|25031|25041|25051|25061|25071|25081|25091|25003|25013|25023|25033|25043|25053|25063|25073|25083|25093", diabetes$ICD9), 1, 0))

table(diabetes$TYPE1DIAB)


#### type 2

diabetes <- diabetes %>%
  mutate(TYPE2DIAB = if_else(
    grepl("E11|O241", ICD10) | 
      grepl("1223|1220", SELFREPORT) |
      grepl("25000|25010|25020|25030|25040|25050|25060|25070|25080|25090|25002|25012|25022|25032|25042|25052|25062|25072|25082|25092", ICD9) |
      grepl("1140868902|1140874646|1140874674|1140874718|1140874744|1140883066|1140884600|1141152590|1141157284|1141168660|1141171646|1141173882|1141189090", MEDICATION) |
      replace_na(HBA1C > 48, FALSE), 
    1, 0
  ))

table(diabetes$TYPE2DIAB)



###### OTHER ILLNESS ###########################################################


colnames(other_illness) <- c("IID", "ICD10", "ICD9", "SELFREPORT", "MEDICATION")


###### arthritis

other_illness <- other_illness %>%
  mutate(ARTHRITIS = if_else(grepl("M05|M06", other_illness$ICD10) | 
                               grepl("1464", other_illness$SELFREPORT)|
                               grepl("714", other_illness$ICD9), 1, 0))


##### afib 

other_illness <- other_illness %>%
  mutate(AFIB = if_else(grepl("I48", other_illness$ICD10) | 
                          grepl("1471|1483", other_illness$SELFREPORT)|
                          grepl("4273|4720", other_illness$ICD9), 1, 0))


##### chronic kidney disease

other_illness <- other_illness %>%
  mutate(KIDNEY_DISEASE = if_else(grepl("N183|N184|N185", other_illness$ICD10) | 
                                    grepl("1192|1519|1609", other_illness$SELFREPORT)|
                                    grepl("5853|5855|5810|5820|5900|V420|V451", other_illness$ICD9), 1, 0))




##### migraine

other_illness <- other_illness %>%
  mutate(MIGRAINE = if_else(grepl("G43|G440|N943", other_illness$ICD10) | 
                              grepl("1265", other_illness$SELFREPORT)|
                              grepl("346", other_illness$ICD9), 1, 0))


##### SLE 

other_illness <- other_illness %>%
  mutate(SLE = if_else(grepl("M32", other_illness$ICD10) | 
                         grepl("1381", other_illness$SELFREPORT)|
                         grepl("7100", other_illness$ICD9), 1, 0))



##### MENTAL ILLNESS

other_illness <- other_illness %>%
  mutate(MENTAL_ILLNESS = if_else(grepl("F03|F068|F09|F20|F22|F23|F259|F28|F29|F31|F39|F53|F333", other_illness$ICD10) | 
                                    grepl("1289|1291", other_illness$SELFREPORT)|
                                    grepl("295|298|296", other_illness$ICD9), 1, 0))



##### ED

other_illness <- other_illness %>%
  mutate(ED = if_else(grepl("N484", other_illness$ICD10) | 
                        grepl("1518", other_illness$SELFREPORT)|
                        grepl("60784", other_illness$ICD9)|
                        grepl("1141168936|1141168948|1141168944|1141168946|1140869100|1140883010", other_illness$MEDICATION), 1, 0))




##################### MEDICATIONS #############################################

medications <- read.csv("TestosteroneCAD/medications_participant.csv")

colnames(medications) <- c("IID", "MEDICATION", "BPMEDS")


medications <- medications %>%
  mutate(HYPERTENSION = if_else(grepl("1140860192|1140860292|1140860696|
                                      1140860728|1140860750|1140860806|1140860882|
                                      1140860904|1140861088|1140861190|1140861276|
                                      1140866072|1140866078|1140866090|1140866102|
                                      1140866108|1140866122|1140866138|1140866156|
                                      1140866162|1140866724|1140866738|1140868618|
                                      1140872568|1140874706|1140874744|1140875808|
                                      1140879758|1140879760|1140879762|1140879802|
                                      1140879806|1140879810|1140879818|1140879822|
                                      1140879826|1140879830|1140879834|1140879842|
                                      1140879866|1140884298|1140888552|1140888556|
                                      1140888560|1140888646|1140909706|1140910442|
                                      1140910614|1140916356|1140923272|1140923336|
                                      1140923404|1140923712|1140926778|1140928226|
                                      1141145660|1141146126|1141152998|1141153026|
                                      1141164276|1141165470|1141166006|1141169516|
                                      1141171336|1141180592|1141180772|1141180778|
                                      1141184722|1141193282|1141194794|1141194810", medications$MEDICATION) | 
                                  grepl("2", medications$BPMEDS), 1, 0))



##### corticosteroids 


medications <- medications %>%
  mutate(CORTICOSTEROIDS = if_else(grepl("1140874790|1140874816|
1140874896.00|1140874930|1140874976|1141145782|1141173346", medications$MEDICATION) , 1, 0))


table(medications$CORTICOSTEROIDS)



##### antipsychotics 


medications <- medications %>%
  mutate(ANTIPSYCHOTICS = if_else(grepl("1140867420|1140867444|1140927956|
                                        1140928916|1141152848|1141153490|
                                        1141169714|1141195974", medications$MEDICATION) , 1, 0))


table(medications$ANTIPSYCHOTICS)







# SOCIODEMOGRAPHIC VARIABLES 

sociodemographics <- read.csv("TestosteroneCAD/sociodemographics_participant.csv")


### smoking


# smoking - this one is more tricky 

sociodemographics$exsmoker <- ifelse(sociodemographics$Ever.smoked...Instance.0 == "1" & sociodemographics$Current.tobacco.smoking...Instance.0=="0" ,"2",NA)
sociodemographics$nonsmoker <- ifelse(sociodemographics$Ever.smoked...Instance.0 == "0" ,"1",NA)
sociodemographics$NUM_CIGS_DAILY <- as.numeric(sociodemographics$Number.of.cigarettes.currently.smoked.daily..current.cigarette.smokers....Instance.0)
sociodemographics$SmokingCategory <- ifelse(sociodemographics$NUM_CIGS_DAILY < 10 & sociodemographics$NUM_CIGS_DAILY > 0, "3",
                                            ifelse(sociodemographics$NUM_CIGS_DAILY >= 10 & sociodemographics$NUM_CIGS_DAILY < 20, "4", 
                                                   ifelse(sociodemographics$NUM_CIGS_DAILY >= 20, "5", NA)))



# Remove NA values and replace with empty strings
sociodemographics$exsmoker[is.na(sociodemographics$exsmoker)] <- ""
sociodemographics$nonsmoker[is.na(sociodemographics$nonsmoker)] <- ""
sociodemographics$SmokingCategory[is.na(sociodemographics$SmokingCategory)] <- ""

# Combine columns into UKBBSMOKING with no white space
sociodemographics$UKBBSMOKING <- paste0(sociodemographics$exsmoker, sociodemographics$nonsmoker, sociodemographics$SmokingCategory)


sociodemographics$UKBBSMOKING <- as.factor(sociodemographics$UKBBSMOKING)

sociodemographics %>%
  mutate(UKBBSMOKING = factor(UKBBSMOKING,
                              levels = c("","1","2","3","4","5")))

sociodemographics <- sociodemographics %>% select(-"nonsmoker", -"SmokingCategory", -"NUM_CIGS_DAILY")





colnames(sociodemographics) <- c("IID", "AGERECRUIT", "MONTHBIRTH", "YEARBIRTH", "DEPRIVATION", "LOST_TO_FOLLOW_UP", "ETHNICITY", 
                                 "BMI", "EVERSMOKED", "SMOKINGSTATUS", "CURRENTSMOKING", "CIGSDAILY", "SBP1", "DATEASSESSMENT", 
                                 "SBP2", "CHOLESTEROL", "HDL", "ILLFATH", "ILLMOTH", "ILLSIBS", "EXSMOKER", "UKBBSMOKING")





##### BLOOD PRESSURE 

invalid_columns <- which(names(sociodemographics) == "" | is.na(names(sociodemographics)))
print(invalid_columns)

names(sociodemographics)[invalid_columns] <- paste0("InvalidName", seq_along(invalid_columns))

sociodemographics <- sociodemographics %>%
  mutate(SBP = (SBP1 + SBP2) / 2)


sociodemographics <- sociodemographics %>%
  rowwise() %>%
  mutate(SBP_SD = sd(c(SBP1, SBP2))) %>%
  ungroup()







##### CHOLESTEROL HDL RATIO


sociodemographics <- sociodemographics %>%
  mutate(CHOLESTEROLTOHDL = CHOLESTEROL / HDL)






##### ethnicity 


sociodemographics <- sociodemographics %>%
  mutate(ETHNICITY_CATEGORY = case_when(
    ETHNICITY %in% c(1, 1001, 1002, 1003, -1, -3) ~ 1,  # White or not stated
    ETHNICITY == 3001 ~ 2,  # Indian
    ETHNICITY == 3002 ~ 3,  # Pakistani
    ETHNICITY == 3003 ~ 4,  # Bangladeshi
    ETHNICITY == 3004 ~ 5,  # Other Asian
    ETHNICITY == 4001 ~ 6,  # Black Caribbean
    ETHNICITY == 4002 ~ 7,  # Black African
    ETHNICITY == 5 ~ 8,  # Chinese
    ETHNICITY == 6 ~ 9,  # Other ethnic group
    TRUE ~ NA_real_  # Default case, if none of the conditions match
  ))



sociodemographics <- sociodemographics %>%
  mutate(ETHNICITY_CATEGORY = factor(ETHNICITY_CATEGORY, levels = 1:9, labels = c(
    "White or not stated",
    "Indian",
    "Pakistani",
    "Bangladeshi",
    "Other Asian",
    "Black Caribbean",
    "Black African",
    "Chinese",
    "Other ethnic group"
  )))




##### ill family member

sociodemographics <- sociodemographics %>%
  mutate(FAMHISTORY = if_else(
    grepl("\\b1\\b", ILLFATH) | grepl("\\b1\\b", ILLMOTH) | grepl("\\b1\\b", ILLSIBS),
    1, 0
  ))








#### selecting the important variables before merging 

medications <- medications %>% select("IID", "HYPERTENSION", "CORTICOSTEROIDS", "ANTIPSYCHOTICS")
diabetes <- diabetes %>% select("IID", "TYPE1DIAB", "TYPE2DIAB")
other_illness <- other_illness %>% select(-c("MEDICATION", "SELFREPORT"))





###### MERGING ALL OF THESE TABLES 

COVARIATES <- merge(diabetes, medications, by = "IID")
COVARIATES <- merge(COVARIATES, other_illness, by = "IID")
COVARIATES <- merge(COVARIATES, sociodemographics, by = "IID")








###### merging with the survival data

COMPLETE_DATA <- merge(surv, COVARIATES, by = "IID")










COMPLETE_DATA <- COMPLETE_DATA %>% select(-c("CIGSDAILY"))
COMPLETE_DATA <- COMPLETE_DATA %>% select(-c("EXSMOKER"))
COMPLETE_DATA <- COMPLETE_DATA %>% select(-c("ICD10", "ICD9"))
COMPLETE_DATA <- COMPLETE_DATA %>% select(-c("LOST_TO_FOLLOW_UP"))
COMPLETE_DATA <- COMPLETE_DATA %>% select(-c("EVERSMOKED", "SMOKINGSTATUS", "CURRENTSMOKING"))
COMPLETE_DATA <- COMPLETE_DATA %>% select(-c("X"))
COMPLETE_DATA <- COMPLETE_DATA %>% select(-c("ILLFATH", "ILLMOTH", "ILLSIBS"))
COMPLETE_DATA <- COMPLETE_DATA %>% select(-c("SBP1", "SBP2"))

any(duplicated(surv_with_dates$IID))
any(duplicated(COMPLETE_DATA$IID))
surv_with_dates <- surv_with_dates[!duplicated(surv_with_dates$IID), ]
COMPLETE_DATA <- COMPLETE_DATA[!duplicated(COMPLETE_DATA$IID), ]



REMOVE_CAD_PRIOR_TO_ASSESSMENT <- merge(surv_with_dates, COMPLETE_DATA, by = "IID")


REMOVE_CAD_PRIOR_TO_ASSESSMENT$DATEASSESSMENT.x <- as.Date(REMOVE_CAD_PRIOR_TO_ASSESSMENT$DATEASSESSMENT.x)

table(REMOVE_CAD_PRIOR_TO_ASSESSMENT$earliest_cad_date_all <= REMOVE_CAD_PRIOR_TO_ASSESSMENT$DATEASSESSMENT.x)

table(REMOVE_CAD_PRIOR_TO_ASSESSMENT$CADBIN.x)



dates <- read.csv("TestosteroneCAD/attendance_participant.csv")

names(dates)[names(dates) == "Participant.ID"] <- "IID"

dat <- merge(REMOVE_CAD_PRIOR_TO_ASSESSMENT, dates, by = "IID")

names(dat)[names(dat) == "Date.of.attending.assessment.centre...Instance.0"] <- "ASSESSMENT_DATE"


cleaned_data <- dat %>%
  filter(is.na(earliest_cad_date_all) | is.na(ASSESSMENT_DATE) | earliest_cad_date_all >= ASSESSMENT_DATE)



cleaned_data <- cleaned_data %>% select(-c("MONTHBIRTH.x"))
cleaned_data <- cleaned_data %>% select(-c("YEARBIRTH.x"))
cleaned_data <- cleaned_data %>% select(-c("DATEASSESSMENT.x"))
cleaned_data <- cleaned_data %>% select(-c("AGERECRUIT.y"))
cleaned_data <- cleaned_data %>% select(-c("LTF"))
cleaned_data <- cleaned_data %>% select(-c("LTFBIN"))
cleaned_data <- cleaned_data %>% select(-c("censdate"))
cleaned_data <- cleaned_data %>% select(-c("cadmonth"))
cleaned_data <- cleaned_data %>% select(-c("cadday"))
cleaned_data <- cleaned_data %>% select(-c("cadyear"))
cleaned_data <- cleaned_data %>% select(-c("censyear"))
cleaned_data <- cleaned_data %>% select(-c("T.y"))
cleaned_data <- cleaned_data %>% select(-c("CADBIN.y"))
cleaned_data <- cleaned_data %>% select(-c("MONTHBIRTH.y"))
cleaned_data <- cleaned_data %>% select(-c("YEARBIRTH.y"))


table(cleaned_data$ASSESSMENT_DATE > cleaned_data$earliest_cad_date_all)




cleaned_data$got_CAD <- !is.na(cleaned_data$earliest_cad_date_all)
cleaned_data$date_comparison <- ifelse(cleaned_data$got_CAD, 
                                       ifelse(cleaned_data$earliest_cad_date_all < cleaned_data$ASSESSMENT_DATE, "Before", "After"),
                                       NA)


# Count the number of participants who never got CAD
num_na_CAD <- sum(is.na(cleaned_data$earliest_cad_date_all))
sum(cleaned_data$CADBIN.x==1)


table(cleaned_data$date_comparison, useNA = "ifany")

write.csv(cleaned_data, file = "COMPLETE_SURV_DATA.csv", row.names = TRUE)



###### THIS SCRIPT ALSO EXISTS IN THE q_risk_time_in_biobank 0907.R script
####### in your local environment 
####### that script contains the data cleaning 
####### in this section of the script we are running the models as this 
####### is the part that cannot be done in R 

install.packages("tidyverse")
install.packages("survival")
install.packages("ggplot2")
install.packages("broom")
install.packages("dplyr")
library(broom)
library(dplyr)
library(tidyverse)
library(survival)
library(ggplot2)


cleaned_data <- read.csv("COMPLETE_SURV_DAT-NEW.csv")


cleaned_data <- cleaned_data %>% select(-c("earliest_cad_date_all"))
cleaned_data <- cleaned_data %>% select(-c("date_comparison"))


complete_data <- cleaned_data[complete.cases(cleaned_data), ]

nrow(complete_data)

qrisksurv <- complete_data

names(qrisksurv)[names(qrisksurv) == "T.x"] <- "T"
names(qrisksurv)[names(qrisksurv) == "CADBIN.x"] <- "CADBIN"


# fixing the ethnicity labelling 



unique(complete_data$ETHNICITY_CATEGORY)

table(complete_data$ETHNICITY_CATEGORY)

complete_data$ETHNICITY_CATEGORY <- as.factor(complete_data$ETHNICITY_CATEGORY)


complete_data$ETHNICITY_CATEGORY <- relevel(complete_data$ETHNICITY_CATEGORY, ref = "White or not stated")



# creating quartiles of the testosterone distribution

quartiles <- quantile(qrisksurv$T, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
qrisksurv$testosterone_quartile <- cut(qrisksurv$T, breaks = quartiles,
                                       labels = c("lower", "lower middle", "upper middle", "upper"),
                                       include.lowest = TRUE)






# age categories 

qrisksurv$age_group <- cut(qrisksurv$AGERECRUIT.x,
                           breaks = c(40, 45, 50, 55, 60, 65, 70, Inf),
                           labels = c("40-45", "45-50", "50-55", "55-60", "60-65", "65-70", "70+"),
                           right = FALSE)


complete_data <- qrisksurv



### smoking variable


complete_data$UKBBSMOKING <- factor(complete_data$UKBBSMOKING, levels = c("1", "2", "3", "4", "5"))
complete_data$UKBBSMOKING <- relevel(complete_data$UKBBSMOKING, ref = "1")

# deficient, sufficient, and high groups 

# Define cutoff points and labels
cut_points <- c(-Inf, 12, 18, 25, 30)
labels <- c("deficient", "sufficient", "high", "very high")

# Create categorical variable for testosterone categories
complete_data$testosterone_category <- cut(complete_data$T, breaks = cut_points, labels = labels, include.lowest = TRUE)

complete_data$testosterone_category <- relevel(complete_data$testosterone_category, ref = "sufficient")




# SURVIVAL ANALYSIS 



# making adjustments - just looking at deficient and sufficient 

complete_data$testosterone_binary <- ifelse(complete_data$T < 12, "deficient", "sufficient")

complete_data$testosterone_binary<- factor(complete_data$testosterone_binary, levels = c("deficient", "sufficient"))

# Reorder levels so that "deficient" is the reference category
complete_data$testosterone_binary <- relevel(complete_data$testosterone_binary, ref = "sufficient")


# Create survival object
CADsurv <- Surv(time = complete_data$timetoEVENT, event = complete_data$CADBIN)

# Fit Kaplan-Meier survival curves
km_fit <- survfit(CADsurv ~ testosterone_binary, data = complete_data)



plot(
  km_fit,
  col = c("blue", "red"),           # Colors for the curves
  lty = 1:2,                        # Line types for the curves
  xlab = "Time to Event",
  ylab = "Survival Probability",
  main = "Kaplan-Meier Survival Curves"
)







# Function to calculate rates and CIs
calculate_rates <- function(complete_data, age_group) {
  group_data <- complete_data %>% filter(age_group == !!age_group)
  person_years <- sum(group_data$timetoEVENT)
  incident_cases <- sum(group_data$CADBIN)
  rate_per_1000 <- (incident_cases / person_years) * 1000
  se_rate <- sqrt(incident_cases) / person_years * 1000
  lower_ci <- rate_per_1000 - 1.96 * se_rate
  upper_ci <- rate_per_1000 + 1.96 * se_rate
  c(Rate_per_1000 = rate_per_1000, Lower_CI = lower_ci, Upper_CI = upper_ci, Person_Years = person_years, Incident_Cases = incident_cases)
}

# List of age groups
age_groups <- levels(complete_data$age_group)

# Calculate for each group
results <- data.frame()
for (age_group in age_groups) {
  rates <- calculate_rates(complete_data, age_group)
  results <- rbind(results, data.frame(Age_Group = age_group, t(rates)))
}

# View results
print(results)

write.csv(results, "cardiovascular_disease_rates.csv", row.names = FALSE)



cox <- coxph(CADsurv~testosterone_binary, data = complete_data)
summary(cox)

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$BMI)
summary(cox) 

cox <- coxph(CADsurv~complete_data$BMI)
summary(cox)

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$TYPE1DIAB)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$TYPE2DIAB)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$HYPERTENSION)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$CORTICOSTEROIDS)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$ANTIPSYCHOTICS)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$ARTHRITIS)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$AFIB)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$KIDNEY_DISEASE)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$MIGRAINE)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$SLE)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$MENTAL_ILLNESS)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$ED)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$AGERECRUIT)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$DEPRIVATION)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$ETHNICITY_CATEGORY)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$BMI)
summary(cox) 


cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$UKBBSMOKING)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$SBP)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$SBP_SD)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$CHOLESTEROLTOHDL)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$FAMHISTORY)
summary(cox) 







# Define a function to extract Cox model results from the summary
extract_cox_summary <- function(model, model_name) {
  summary_model <- summary(model)
  coef <- summary_model$coefficients
  confint <- summary_model$conf.int
  results <- data.frame(
    model = model_name,
    hazard_ratio = coef[, "exp(coef)"],
    p.value = coef[, "Pr(>|z|)"],
    conf.low = confint[, "lower .95"],
    conf.high = confint[, "upper .95"]
  )
  return(results)
}

# Define a list of models with their names
models <- list(
  "testosterone_binary" = coxph(CADsurv ~ testosterone_binary, data = complete_data),
  "testosterone_binary + BMI" = coxph(CADsurv ~ testosterone_binary + BMI, data = complete_data),
  "testosterone_binary + TYPE1DIAB" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$TYPE1DIAB),
  "testosterone_binary + TYPE2DIAB" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$TYPE2DIAB),
  "testosterone_binary + HYPERTENSION" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$HYPERTENSION),
  "testosterone_binary + CORTICOSTEROIDS" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$CORTICOSTEROIDS),
  "testosterone_binary + ANTIPSYCHOTICS" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$ANTIPSYCHOTICS),
  "testosterone_binary + ARTHRITIS" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$ARTHRITIS),
  "testosterone_binary + AFIB" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$AFIB),
  "testosterone_binary + KIDNEY_DISEASE" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$KIDNEY_DISEASE),
  "testosterone_binary + MIGRAINE" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$MIGRAINE),
  "testosterone_binary + SLE" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$SLE),
  "testosterone_binary + MENTAL_ILLNESS" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$MENTAL_ILLNESS),
  "testosterone_binary + ED" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$ED),
  "testosterone_binary + AGERECRUIT" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$AGERECRUIT),
  "testosterone_binary + DEPRIVATION" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$DEPRIVATION),
  "testosterone_binary + ETHNICITY_CATEGORY" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$ETHNICITY_CATEGORY),
  "testosterone_binary + UKBBSMOKING" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$UKBBSMOKING),
  "testosterone_binary + SBP" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$SBP),
  "testosterone_binary + SBP_SD" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$SBP_SD),
  "testosterone_binary + CHOLESTEROLTOHDL" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$CHOLESTEROLTOHDL),
  "testosterone_binary + FAMHISTORY" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$FAMHISTORY)
)

# Extract results for each model and combine them into a single table
combined_results <- bind_rows(lapply(names(models), function(model_name) {
  model <- models[[model_name]]
  extract_cox_summary(model, model_name)
}))

# Print the combined results
print(combined_results)

# Save the results to a CSV file (optional)
write.csv(combined_results, file = "separate_mediators.csv", row.names = TRUE)













##### ADJUSTING FOR EVERYTHING 

# Fit Cox proportional hazards model with multiple predictors
cox <- coxph(
  formula = CADsurv ~ testosterone_binary + BMI + TYPE1DIAB + TYPE2DIAB + HYPERTENSION +
    CORTICOSTEROIDS + ANTIPSYCHOTICS + ARTHRITIS + AFIB + KIDNEY_DISEASE +
    MIGRAINE + SLE + MENTAL_ILLNESS + ED + AGERECRUIT.x + DEPRIVATION +
    ETHNICITY_CATEGORY + UKBBSMOKING + SBP + SBP_SD + CHOLESTEROLTOHDL +
    FAMHISTORY,
  data = complete_data
)

# Summarize the model
summary(cox)



cox_summary <- summary(cox)


results <- data.frame(
  exp_coef = cox_summary$coefficients[, "exp(coef)"],
  p_value = cox_summary$coefficients[, "Pr(>|z|)"],
  lower_ci = cox_summary$conf.int[, 1],  # Lower 95% CI
  upper_ci = cox_summary$conf.int[, 2]   # Upper 95% CI
)

write.csv(results, file = "cox_model_results.csv", row.names = TRUE)










### now having a LOOK AT THE DISTRIBUTION BY DECILES ############################

complete_data <- complete_data %>%
  mutate(testosterone_decile = ntile(T, 10))

complete_data <- complete_data %>%
  mutate(testosterone_decile = ntile(T, 10)) %>%
  mutate(testosterone_decile = factor(testosterone_decile)) %>%
  mutate(testosterone_decile = fct_relevel(testosterone_decile, "7"))


# Check the factor levels to confirm releveling
levels(complete_data$testosterone_decile)


decile_ranges <- complete_data %>%
  group_by(testosterone_decile) %>%
  summarise(
    min_value = min(T),
    max_value = max(T),
    .groups = 'drop'
  )


# Create survival object
CADsurv <- Surv(time = complete_data$timetoEVENT, event = complete_data$CADBIN)

# Fit Kaplan-Meier survival curves
km_fit <- survfit(CADsurv ~ testosterone_decile + AGERECRUIT.x, data = complete_data)

print(km_fit)



cox_model <- coxph(CADsurv~testosterone_decile+AGERECRUIT.x, data = complete_data)
summary_cox <- summary(cox_model)




# Extracting the coefficients table
coef_table <- summary_cox$coefficients

# Extracting the confidence intervals
conf_int <- summary_cox$conf.int

# Creating the dataframe
results_df <- data.frame(
  coef = coef_table[, "coef"],
  exp_coef = coef_table[, "exp(coef)"],
  se_coef = coef_table[, "se(coef)"],
  z = coef_table[, "z"],
  Pr_z = coef_table[, "Pr(>|z|)"],
  exp_coef_lower_95 = conf_int[, "lower .95"],
  exp_coef_upper_95 = conf_int[, "upper .95"]
)

# Adding row names
rownames(results_df) <- rownames(coef_table)

# Display the dataframe
print(results_df)



write.csv(results_df, file = "decile_stratified.csv", row.names = TRUE)


















### now having a LOOK AT THE DISTRIBUTION BY QUARTILES ############################


complete_data$testosterone_quartile <- factor(complete_data$testosterone_quartile, 
                                              levels = c("upper", "lower", "lower middle", "upper middle"))

# Check the factor levels to confirm releveling
levels(complete_data$testosterone_quartile)


# Create survival object
CADsurv <- Surv(time = complete_data$timetoEVENT, event = complete_data$CADBIN)

# Fit Kaplan-Meier survival curves
km_fit <- survfit(CADsurv ~ testosterone_quartile, data = complete_data)





cox_model <- coxph(CADsurv~testosterone_quartile + AGERECRUIT.x, data = complete_data)
summary_cox <- summary(cox_model)




# Extracting the coefficients table
coef_table <- summary_cox$coefficients

# Extracting the confidence intervals
conf_int <- summary_cox$conf.int

# Creating the dataframe
results_df <- data.frame(
  coef = coef_table[, "coef"],
  exp_coef = coef_table[, "exp(coef)"],
  se_coef = coef_table[, "se(coef)"],
  z = coef_table[, "z"],
  Pr_z = coef_table[, "Pr(>|z|)"],
  exp_coef_lower_95 = conf_int[, "lower .95"],
  exp_coef_upper_95 = conf_int[, "upper .95"]
)

# Adding row names
rownames(results_df) <- rownames(coef_table)

# Display the dataframe
print(results_df)



write.csv(results_df, file = "quartile_stratified.csv", row.names = TRUE)





####### READING IN THE RESULTS FROM THE RAP TO CREATE SOME PLOTS ######

person_years <- read.csv("TestosteroneCAD/Data_from_RAP/cardiovascular_disease_rates.csv")
cox_model_all_mediators <- read.csv("TestosteroneCAD/Data_from_RAP/cox_model_results.csv")
cox_model_stratified <- read.csv("TestosteroneCAD/Data_from_RAP/combined_cox_model_results.csv")
quartile_stratified <- read.csv("TestosteroneCAD/Data_from_RAP/quartile_stratified.csv")
clinical <- read.csv("TestosteroneCAD/Data_from_RAP/clinical_T_cox.csv")
decile_stratified <- read.csv("TestosteroneCAD/Data_from_RAP/decile_stratified.csv")

##### PLOTTING THE CAD HAZARDS FOR TESTOSTERONE QUARTILES IN THE WHOLE 
#### POPULATION 

names(quartile_stratified)[names(quartile_stratified) == "X"] <- "T_QUARTILE"
print(quartile_stratified)


# Reorder the factor levels for T_QUARTILE
quartile_stratified$T_QUARTILE <- factor(quartile_stratified$T_QUARTILE,
                                         levels = c("lower", "lower middle", "upper middle", "upper"))

# Custom colors for testosterone quartiles
testosterone_colors <- c("lower" = "#66c2a5", "lower middle" = "#fc8d62", "upper middle" = "#8da0cb", "upper" = "#e78ac3")

# Create the forest plot without a legend
plot <- ggplot(quartile_stratified, aes(x = T_QUARTILE, y = exp_coef, ymin = exp_coef_lower_95, ymax = exp_coef_upper_95, color = T_QUARTILE)) +
  geom_point(position = position_dodge(width = 0.3), size = 4) +
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.2) +
  labs(title = "Hazard Ratios of CAD by Testosterone Quartile",
       x = "Testosterone Quartile",
       y = "Hazard Ratio of CAD") +
  scale_color_manual(values = testosterone_colors, guide = FALSE) +  # Use custom colors without legend
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size = 14)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  # Add a line at HR = 1 (null effect)
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))


print(plot)

ggsave("forest_plot_quartiles.png", plot = plot, width = 8, height = 6, dpi = 300)





##### PLOTTING THE CAD HAZARDS FOR TESTOSTERONE DECILES IN THE WHOLE 
#### POPULATION 
head(decile_stratified)
print(decile_stratified)

decile_stratified$X <- factor(decile_stratified$X,
                                         levels = c("testosterone_decile1", "testosterone_decile2", "testosterone_decile3", "testosterone_decile4", "testosterone_decile5",
                                                    "testosterone_decile6", "testosterone_decile7", "testosterone_decile8", "testosterone_decile9", "testosterone_decile10"))

print(decile_stratified)


# Create the forest plot without a legend
plot <- ggplot(decile_stratified, aes(x = X, y = exp_coef, ymin = exp_coef_lower_95, ymax = exp_coef_upper_95)) +
  geom_point(position = position_dodge(width = 0.3), size = 4) +
  geom_errorbar(position = position_dodge(width = 0.6), width = 0.2) +
  labs(title = "Hazard Ratios of CAD by Testosterone Decile",
       x = "Testosterone Decile",
       y = "Hazard Ratio of CAD") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size = 14)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  # Add a line at HR = 1 (null effect)
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))


print(plot)

ggsave("forest_plot_quartiles.png", plot = plot, width = 8, height = 6, dpi = 300)






########### CREATING A FOREST PLOT BASED ON CLINICAL VALUES OF TESTOSTERONE 
########## AS OPPOSED TO THE QUARTILES IN THE DISTRIBUTION 

 
print(clinical)


names(clinical)[names(clinical) == "X"] <- "clinical_T"


# Ensure clinical_T is a factor with the correct order of levels
clinical$clinical_T <- factor(clinical$clinical_T, levels = c("very-low", "low", "standard", "high"))

# Custom colors for clinical T categories in desired order
clinical_colors <- c("very-low" = "#fc8d62", "low" = "#66c2a5", "standard" = "#8da0cb", "high" = "#e78ac3")

# Create the forest plot with reordered levels
plot <- ggplot(clinical, aes(x = clinical_T, y = exp_coef, ymin = exp_coef_lower_95, ymax = exp_coef_upper_95, color = clinical_T)) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_errorbar(position = position_dodge(width = 0.3), width = 0.2) +
  labs(title = "Hazard Ratios of Incident CAD by Clinical Testosterone Levels",
       x = "Clinical Testosterone Level",
       y = "Hazard Ratio",
       color = "Clinical Testosterone Level") +
  scale_color_manual(values = clinical_colors) +  # Use custom colors with reordered levels
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  # Add a line at HR = 1 (null effect)
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),  # Remove horizontal gridlines for cleaner look
        axis.text.x = element_text(angle = 45, hjust = 1)) 

# Print the plot
print(plot)



